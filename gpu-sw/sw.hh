#include <cuda_runtime.h>
//
#include "core.hh"
#include "seqs.hh"
#include "cuda-utils.hh"
#include<iostream>

/*! sw alg.
 *
 * usage:
 *
 * 1.  use eccl::seqs_builder to build a seq batch (in host memory).
 *     call .push_back(seq) to add a sequence.
 *     construct eccl::seqs after sequences are added.
 * 2.  use eccl::sw_buffers to manage device buffers.
 *     call .init(band_size) once to init.
 *     call [res, cigar_buf, cigar_s]=.run(query, target) to run.
 *     call eccl::format_sw_result() to format output to a buffer.
 *     eccl::sw_buffers can be reused.
 *
 * see check_xna() in sw-test.cc for simple examples.
 *
 * NOTES
 * - correct compared to gasal2.
 * - just use the alternative API (eccl::sw_buffers).
 *   not removing existing async stuff.
 * - why no batch?
 *   input is not a single seq, not all the seqs, then it *is* a batch.
 * - pass nullptr if no stream is needed.
 * - the sw part may need further work, i'll do code clean up after that.
 *   print and time related calls are silenced.
 * - non-relevant (and relevant) codes are in a separate dir.
 *   just leave them there and use the API.
 *
 *
 * some rough benchmarks:
 *
 * gasal2 no output: total 0m4.430s,  kernel 1.03s
 * gasal2:           total 0m12.239s, kernel 8.77s
 * eccl:             total 0m4.289s,  kernel 2.4s
 * eccl banded(16):  total 0m3.632s,  kernel 1.8s
 *
 * eccl tests are run asynchronously.
 * there are noticeable overheads in data copying, i'll improve that later.
 *
 */

namespace eccl {

struct sw_result {
	int score;
	uint2 query_se;
	uint2 target_se;
};

/*! global/fixed stuff */
struct sw_gpu_opts {
	unsigned int num_blocks;
	unsigned int num_threads;
	unsigned int band_size;
	dev_ptr<signed short> subst_mat;
};

/*! per-batch stuff */
struct sw_gpu_args {
	cudaStream_t stream;
	unsigned int work_size;
	unsigned int max_query_len;
	dev_ptr<unsigned int> target_buf;
	dev_ptr<ulong2> target_se;
	dev_ptr<unsigned int> query_buf;
	dev_ptr<ulong2> query_se;
	dev_ptr<int2> tmp_scores;
	dev_ptr<unsigned int> tmp_traceback_buf;
	dev_ptr<unsigned long> tmp_traceback_s;
	dev_ptr<sw_result> results;
	dev_ptr<unsigned int> cigar_buf;
	dev_ptr<unsigned long> cigar_s;
};

inline unsigned int estimate_cigar_size(unsigned int query_len, unsigned int target_len) {
	// XXX upper bound???
	//return query_len/4;
	//return 48;
	//return 20;
	return 48;
	// (450)/16 == 28
	// no_cigar(flat_traceback): 1.31s
	// 8: 1.48s
	// 16: 1.45s
	// 32: 1.50s
	// 64: 1.78s
	// 128: 3.0s
}

template<seq_type Type>
void sw_gpu(const sw_gpu_opts& opts, sw_gpu_args&& args);

template<seq_type Type>
struct sw_buffers {
	unsigned int band_size;
	eccl::dev_vector<signed short> subst_mat;
	eccl::dev_vector<unsigned int> query_buf;
	eccl::dev_vector<ulong2> query_se;
	eccl::dev_vector<unsigned int> target_buf;
	eccl::dev_vector<ulong2> target_se;
	eccl::dev_vector<int2> tmp_scores;
	eccl::dev_vector<ulong> tmp_traceback_s;
	eccl::dev_vector<unsigned int> tmp_traceback_buf;
	eccl::dev_vector<eccl::sw_result> results;
	eccl::dev_vector<ulong> cigar_s;
	eccl::dev_vector<unsigned int> cigar_buf;

	/*! results.  use after run */
	eccl::host_vector<eccl::sw_result> results_host;
	eccl::host_vector<unsigned int> cigar_buf_host;
	eccl::host_vector<ulong> cigar_s_host;

	template<typename ScoreOpts>
	void init(unsigned int band_size, const ScoreOpts& scores) {
		this->band_size=band_size;
		eccl::host_vector<signed short> subst_mat_host;
		eccl::fill_score_table<Type>(subst_mat_host, scores);
		subst_mat.malloc_h2d(subst_mat_host, nullptr);
	}
	void init(unsigned int band_size) {
		this->band_size=band_size;
		eccl::host_vector<signed short> subst_mat_host;
		eccl::fill_score_table<Type>(subst_mat_host);
		subst_mat.malloc_h2d(subst_mat_host, nullptr);
	}
	inline void run(const seqs<Type>& query, const seqs<Type>& target);
};

template<seq_type Type>
inline void sw_buffers<Type>::run(const seqs<Type>& query, const seqs<Type>& target) {
	sw_gpu_opts opts;
	opts.num_blocks=40;
	opts.num_threads=128;
	opts.band_size=this->band_size;
	opts.subst_mat=this->subst_mat;

	sw_gpu_args args;
	args.stream=nullptr;
	assert(query.size()==target.size());
	args.work_size=query.size();
	args.max_query_len=eccl::padded_len<Type>(query.max_seq_len());

	eccl::host_vector<unsigned int> query_buf_host;
	query_buf_host.insert(query_buf_host.end(), query.buf(), query.buf()+query.buf_size());
	query_buf.try_malloc_h2d(query_buf_host, nullptr);
	eccl::host_vector<ulong2> query_se_host;
	for(unsigned int i=0; i<query.size(); ++i) {
		auto begin=query.begin(i);
		auto seq=query[i];
		query_se_host.push_back({begin+seq.begin, begin+seq.end});
	}
	query_se.try_malloc_h2d(query_se_host, nullptr);

	eccl::host_vector<unsigned int> target_buf_host;
	target_buf_host.insert(target_buf_host.end(), target.buf(), target.buf()+target.buf_size());
	target_buf.try_malloc_h2d(target_buf_host, nullptr);
	eccl::host_vector<ulong2> target_se_host;
	for(unsigned int i=0; i<query.size(); ++i) {
		auto begin=target.begin(i);
		auto seq=target[i];
		target_se_host.push_back({begin+seq.begin, begin+seq.end});
	}
	target_se.try_malloc_h2d(target_se_host, nullptr);

	tmp_scores.try_malloc(opts.num_threads*opts.num_blocks*args.max_query_len);
	eccl::host_vector<ulong> tmp_traceback_s_host;
	tmp_traceback_s_host.push_back(0);
	for(unsigned int i=0; i<opts.num_threads*opts.num_blocks; ++i) {
		std::size_t n{0};
		for(unsigned int j=i; j<query.size(); j+=opts.num_threads*opts.num_blocks) {
			auto s=target[j].len();
			n=std::max<std::size_t>(n, s);
		}
		n=args.max_query_len*((n+31)/32)*4;
		tmp_traceback_s_host.push_back(tmp_traceback_s_host.back()+n);
	}
	tmp_traceback_s.try_malloc_h2d(tmp_traceback_s_host, nullptr);
	tmp_traceback_buf.try_malloc(tmp_traceback_s_host.back());

	results.try_malloc(query.size());

	cigar_s_host.clear();
	cigar_s_host.push_back(0);
	for(unsigned int i=0; i<query.size(); ++i) {
		auto len=eccl::estimate_cigar_size(query[i].len(),
				target[i].len());
		cigar_s_host.push_back(cigar_s_host.back()+(len+3)/4+1);
	}
	cigar_s.try_malloc_h2d(cigar_s_host, nullptr);
	cigar_buf.try_malloc(cigar_s_host.back());

	args.target_buf=target_buf;
	args.target_se=target_se;
	args.query_buf=query_buf;
	args.query_se=query_se;
	args.tmp_scores=tmp_scores;
	args.tmp_traceback_buf=tmp_traceback_buf;
	args.tmp_traceback_s=tmp_traceback_s;
	args.results=results;
	args.cigar_buf=cigar_buf;
	args.cigar_s=cigar_s;
	sw_gpu<Type>(opts, std::move(args));
	cudaGetLastError(), eccl::check_cuda{"align"};
	cigar_buf_host.resize(cigar_buf.size());
	cigar_buf.d2h(cigar_buf_host, nullptr);
	results_host.resize(results.size());
	results.d2h(results_host, nullptr);
	cudaStreamSynchronize(nullptr), check_cuda{"sync stream"};
};

std::string_view format_sw_result(std::array<char, 4096>& buf,
		const sw_result& align,
		std::string_view name_query,
		std::string_view name_target,
		const unsigned int* cigar,
		unsigned int cigar_len);

}
