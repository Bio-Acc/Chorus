
#include "sw.hh"
#include "sw-lib.hh"

#include <iostream>
#include <charconv>
#include <typeinfo>
using namespace glf;
// struct glf::sw_handle {
	// eccl::sw_buffers<eccl::seq_type::xna> buffs;
// };
template<eccl::seq_type T>
struct glf::sw_handle {
	eccl::sw_buffers<T> buffs;
};


template<eccl::seq_type T>
sw_handle<T>* glf::sw_create(const sw_opts& opts) {
	auto res=std::make_unique<sw_handle<T>>();
	res->buffs.init(0, opts.scores);
	return res.release();
}

template<eccl::seq_type T>
sw_handle<T>* glf::sw_create() {
	auto res=std::make_unique<sw_handle<T>>();
	res->buffs.init(0);
	return res.release();
}

template<eccl::seq_type T>
void glf::sw_destroy(sw_handle<T>* hdl) {
	delete hdl;
}


template<eccl::seq_type T>
std::vector<sw_align> glf::sw_batch(sw_handle<T>* hdl,
		const std::vector<std::string_view>& query,
		const std::vector<std::string_view>& target,
		const sw_batch_opts& opts) {
	auto gen_seq_enc=[](auto& raw) {
		eccl::seqs_builder<T> builder{0, 0, 0};
		for(auto s: raw)
			builder.template push_back<false>(s);
		return eccl::seqs<T>{std::move(builder)};
	};
	auto query_enc=gen_seq_enc(query);
	auto target_enc=gen_seq_enc(target);
	auto& buffs=hdl->buffs;
	buffs.run(query_enc, target_enc);
	
	std::vector<sw_align> results;
	results.reserve(query.size());
	for(unsigned int i=0; i<query.size(); ++i) {
		auto& res=results.emplace_back();
		auto& res0=buffs.results_host[i];
		auto off=buffs.cigar_s_host[i];
		res.sw_score=res0.score;
		res.ref_begin=res0.target_se.x;
		res.ref_end=res0.target_se.y+1;
		res.query_begin=res0.query_se.x;
		res.query_end=res0.query_se.y+1;
		auto cigar=&buffs.cigar_buf_host[off+1];
		auto cigar_len=buffs.cigar_buf_host[off];

		res.cigar_len=0;
		auto put_cigar=[&res](unsigned int cnt, unsigned int op) {
			if(res.cigar_len<MAX_NUM_CIGAR) {
				res.cigar[res.cigar_len]=(cnt<<2)|op;
			}
			res.cigar_len+=1;
		};

		unsigned int prev_cnt{0};
		unsigned char prev_op{'\x00'};
		auto nn=cigar_len;
		unsigned int mmcnt=0;
		while(nn>0) {
			--nn;
			unsigned int v=(cigar[nn/4]>>((nn%4)*8))&0xff;
			auto op="\x00\x03\x02\x00"[v&0b11];
			mmcnt+="\x01\x01\x01\x00"[v&0b11];
			v=(v>>2)+1;
			if(op==prev_op) {
				prev_cnt+=v;
			} else {
				if(prev_cnt>0) {
					put_cigar(prev_cnt, prev_op);
				}
				prev_cnt=v;
				prev_op=op;
			}
		}
		if(prev_cnt>0) {
			put_cigar(prev_cnt, prev_op);
		}
		if(res.cigar_len>MAX_NUM_CIGAR)
			res.sw_score=0;
		res.mismatch=mmcnt;
		if(false) {
			std::array<char, 4096> buf;
			std::cerr<<"seq1 "<<query[i]<<" "<<query[i].size();
			std::cerr<<"\nseq2 "<<target[i]<<" "<<target[i].size();
			auto r=format_sw_result(buf, res0, "q", "t", cigar, cigar_len);
			std::cerr<<"\n"<<r<<"\n";
		}
	}
	return results;
}



// static void cigar_to_string(const sw_align& alignment, std::string_view& query, std::string_view& ref, std::vector<int>& q_res, std::vector<size_t>& r_res) {
//   q_res.clear();
//   r_res.clear();
//   int qp = alignment.query_begin;
//   int rp = alignment.ref_begin;
//   for (int c=0;c<alignment.cigar_len;c++)
//   {
//     int cigar_int = cigar_int_to_len(alignment.cigar[c]);
//     char cigar_op = cigar_int_to_op(alignment.cigar[c]);
//     // cout << cigar_int << cigar_op<<endl;
//     if (cigar_op=='M'||cigar_op=='X'){
//       // q_res.append(query+qp, query+qp+cigar_int);
//       for (int i=0;i<cigar_int;i++)
//       {
//         q_res.push_back(qp+i);
//         r_res.push_back(rp+i);
//       }
//       // r_res.append(ref+rp, ref+rp+cigar_int);
//       qp+=cigar_int;
//       rp+=cigar_int;
//     }
//     else if (cigar_op=='I')
//     {
//       // q_res.append(query+qp, query+qp+cigar_int);
//       for (int i=0;i<cigar_int;i++)
//       {
//         q_res.push_back(qp+i);
//         r_res.push_back(-1);
//       }
//       qp+=cigar_int;
//     }
//     else if (cigar_op=='D')
//     {
//       // r_res.append(ref+rp, ref+rp+cigar_int);
//       for (int i=0;i<cigar_int;i++)
//       {
//         r_res.push_back(rp+i);
//         q_res.push_back(-1);
//       }
//       rp+=cigar_int;
//     }
//     else if (cigar_op=='S')
//     {
//       // qp+=cigar_int;
//     }
//     else{
//       std::cout <<"error cigar:"<<cigar_op<<std::endl;
//       exit(-1);
//     }
//   }
// }

template<typename T>
static void inspect(const eccl::dev_vector<T>& obj, sw_stats& stats) {
	stats.dev_mem_usage+=obj.size()*sizeof(T);
}
template<typename T>
static void inspect(const eccl::host_vector<T>& obj, sw_stats& stats) {
	stats.host_mem_usage+=obj.size()*sizeof(T);
}
template<eccl::seq_type Type>
static void inspect(const eccl::sw_buffers<Type>& obj, sw_stats& stats) {
	inspect(obj.subst_mat, stats);
	inspect(obj.query_buf, stats);
	inspect(obj.query_se, stats);
	inspect(obj.target_buf, stats);
	inspect(obj.target_se, stats);
	inspect(obj.tmp_scores, stats);
	inspect(obj.tmp_traceback_s, stats);
	inspect(obj.tmp_traceback_buf, stats);
	inspect(obj.results, stats);
	inspect(obj.cigar_s, stats);
	inspect(obj.cigar_buf, stats);
	inspect(obj.results_host, stats);
	inspect(obj.cigar_buf_host, stats);
	inspect(obj.cigar_s_host, stats);
}

template<eccl::seq_type T>
static void inspect(const sw_handle<T>& obj, sw_stats& stats) {
	inspect(obj.buffs, stats);
}

template<eccl::seq_type T>
sw_stats sw_inspect(sw_handle<T>* hdl) {
	assert(hdl);
	sw_stats stats;
	stats.dev_mem_usage=0;
	stats.host_mem_usage=0;
	inspect(*hdl, stats);
	return stats;
}

inline static char* format_impl(char* p, char* q, const char* n) {
	while(*n) {
		*p=*n;
		++p;
		++n;
	}
	return p;
}
inline static char* format_impl(char* p, char* q, char n) {
	*p=n;
	++p;
	return p;
}
inline static char* format_impl(char* p, char* q, std::string_view n) {
	for(unsigned int i=0; i<n.size(); ++i) {
		*p=n[i];
		++p;
	}
	return p;
}
template<typename T, typename=std::enable_if_t<std::is_integral<T>::value>>
inline static char* format_impl(char* p, char* q, T n) {
	auto res=std::to_chars(p, q, n);
	if(res.ec!=std::errc{})
		throw;
	p=res.ptr;
	return p;
}

std::string_view eccl::format_sw_result(std::array<char, 4096>& buf,
		const eccl::sw_result& align,
		std::string_view name_query,
		std::string_view name_target,
		const unsigned int* cigar,
		unsigned int cigar_len) {
	auto put=[&buf](char* p, auto n) ->char* {
		return format_impl(p, &buf[4096], n);
	};
	auto p=&buf[0];
	p=put(p, "query_name=");
	p=put(p, name_query);
	p=put(p, "\ttarget_name=");
	p=put(p, name_target);
	p=put(p, "\tscore=");
	p=put(p, align.score);
	p=put(p, "\tquery_batch_start=");
	p=put(p, align.query_se.x);
	p=put(p, "\ttarget_batch_start=");
	p=put(p, align.target_se.x);
	p=put(p, "\tquery_batch_end=");
	p=put(p, align.query_se.y);
	p=put(p, "\ttarget_batch_end=");
	p=put(p, align.target_se.y);
	p=put(p, "\tCIGAR=");
	// cigar is slow
	//continue;

	//fprintf(stderr, "alignment: %zd %u\n", oo, nn);
	std::size_t prev_cnt{0};
	char prev_op{'\x00'};
	auto nn=cigar_len;
	while(nn>0) {
		--nn;
		unsigned int v=(cigar[nn/4]>>((nn%4)*8))&0xff;
		auto op="XIDM"[v&0b11];
		v=(v>>2)+1;
		if(op==prev_op) {
			prev_cnt+=v;
		} else {
			if(prev_cnt>0) {
				p=put(p, prev_cnt);
				p=put(p, prev_op);
			}
			prev_cnt=v;
			prev_op=op;
		}
	}
	if(prev_cnt>0) {
		p=put(p, prev_cnt);
		p=put(p, prev_op);
	}
	p=put(p, '\n');
	std::size_t l=p-&buf[0];
	return {&buf[0], l};
}
