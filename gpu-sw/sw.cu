#include <cassert>
#include <chrono>
#include <vector>
#include <charconv>

#include <cuda_runtime.h>

#include "cuda-utils.hh"
#include "sw.hh"
#include "core.hh"



#include <iostream>

using check_cuda=eccl::check_cuda;

template<eccl::seq_type Type>
__global__ void align_kernel(
		/* Req/Reuse/Share/Combine Input Output index */
		/* *!!! * ! */const unsigned int work_size,
		/* *!!! * ! */const unsigned int max_query_len,
		const unsigned int band_size,
		/* ***! * ! */const signed short* subst_mat,
		/* *..! * ! */const unsigned int* target_buf,
		/* *..! * ! work */const ulong2* target_se,
		/* *!!! * ! */const unsigned int* query_buf,
		/* *!!! * ! work */const ulong2* query_se,
		/* *!!! ! ! */int2* tmp_scores,
		/* *!!! ! ! */unsigned int* tmp_traceback_buf,
		/* *!!! * ! thr */const ulong* tmp_traceback_s,
		/* *!!! ! * */eccl::sw_result* results,
		/* *!!! ! * */unsigned int* cigar_buf,
		/* *!!! * ! work */const unsigned long* cigar_s,
		int) {

	const unsigned int thr_id=threadIdx.x+blockIdx.x*blockDim.x;
	const unsigned int debug=0;
	//const unsigned int thr_num=blockDim.x*gridDim.x;
	constexpr uint2 tile_size={8, 8*2};
	const int gap_ext=subst_mat[(Type==eccl::seq_type::xna?8*8:32*32)+1];
	const int gap_open=subst_mat[(Type==eccl::seq_type::xna?8*8:32*32)+0]-gap_ext;

	if(thr_id>=work_size)
		return;
	const auto work_id=thr_id;

	auto max_prev_row=&tmp_scores[thr_id*max_query_len];
	auto t_offset=tmp_traceback_s[thr_id];
	[[maybe_unused]] unsigned int t_size=tmp_traceback_s[thr_id+1]-t_offset;
	assert(t_offset%4==0);
	auto mat2=&tmp_traceback_buf[t_offset];

	struct {
		unsigned int ii;
		unsigned int jj;
		int score;
		unsigned int i0, j0;
	} _top;
	auto a_offset=query_se[work_id].x;
	assert(a_offset<query_se[work_id].y);
	unsigned int a_size=query_se[work_id].y-a_offset;
	auto b_offset=target_se[work_id].x;
	assert(b_offset<target_se[work_id].y);
	unsigned int b_size=target_se[work_id].y-b_offset;
	assert(max_query_len*b_size/8<=t_size);
	if(debug)
	printf("@%u sizes: %u %u\n", thr_id, a_size, b_size);

	_top.score=-9999;

	const unsigned int half_band_width=band_size?band_size:1024*1024*1024;
	auto get_query_range=[a_size,b_size,half_band_width](unsigned int j0, unsigned int j1) ->uint2 {
		unsigned int lb=0, ub=a_size;
		if(a_size<b_size) {
			// 0: [0 hbw)
			// 1: [0 hbw+1)
			// ...
			// b_size-2: [a_size-hbw-1 a_size)
			// b_size-1: [a_size-hbw a_size)
			if(ub>half_band_width+j1)
				ub=half_band_width+j1;
			if(a_size+j0+1>b_size+half_band_width)
				lb=a_size+j0+1-b_size-half_band_width;
		} else {
			// similarly
			if(j0+1>half_band_width)
				lb=j0-half_band_width+1;
			if(half_band_width+j1<b_size)
				ub=a_size-b_size+half_band_width+j1;
		}
		return {lb, ub};
	};
	auto decode_conds=[](unsigned int ec) {
		return (ec*11)%32;
	};
	auto encode_conds=[decode_conds](unsigned char c1, unsigned char c2, unsigned char c3) ->unsigned int {
		/* each c* can be {0, 1, 2};  thus 27 possible cases;  but only 8 casses can occur.
		   encode the 8 cases with 4 bits to save space.
		 */
		auto c=(c1*3u+c2)*3u+c3;
		auto ec=(c*3)%32;
		assert(ec<16);
		assert(decode_conds(ec)==c);
		return ec;
	};

	for(unsigned int i=0; i<max_query_len; ++i) {
		max_prev_row[i]=make_int2(0, gap_open);
	}
	for(unsigned int j=0; j<b_size; j+=tile_size.y) {
		int2 max_prev_col[tile_size.y];
		eccl::code<Type> b_cache[tile_size.y];
		for(unsigned int jj=0; jj<tile_size.y; ++jj) {
			b_cache[jj]=(j+jj<b_size)?eccl::get_code<Type>(target_buf,
					b_offset+j+jj):eccl::code<Type>{};
			max_prev_col[jj]=make_int2(0, gap_open);
		}
		auto qrange=get_query_range(j, j+tile_size.y-1);
		qrange.x=qrange.x/tile_size.x*tile_size.x;
		if(debug)
		printf("@%u i range: %u %u %u\n", thr_id, j, qrange.x, qrange.y);
		for(unsigned int i=qrange.x; i<qrange.y; i+=tile_size.x) {
			unsigned int tb_batch[tile_size.y/8*tile_size.x];
			eccl::code<Type> a_cache[tile_size.x];
			for(unsigned int ii=0; ii<tile_size.x; ++ii) {
				a_cache[ii]=(i+ii<a_size)?eccl::get_code<Type>(query_buf,
						a_offset+i+ii):eccl::code<Type>{};
			}
// per-tile
for(unsigned int ii=0; ii<tile_size.x; ++ii) {
	auto max_prev_row_i_ii=max_prev_row[i+ii];
	for(unsigned int jj=0; jj<tile_size.y; ++jj) {
		int3 score;
		eccl::pair<Type> bp{a_cache[ii], b_cache[jj]};
		auto bp_score=subst_mat[bp.value];
		score.x=max(max_prev_row_i_ii.x+bp_score, 0);
		score.y=max_prev_row_i_ii.y+gap_ext;
		score.z=max_prev_col[jj].y+gap_ext;

		max_prev_row_i_ii.x=max_prev_col[jj].x;
		/* comparison conditions for backtrack */
		unsigned char cond1=0, cond2=0, cond3=0;
		max_prev_col[jj].x=score.x;
		if(max_prev_col[jj].x<score.y) {
			max_prev_col[jj].x=score.y;
			cond1=1;
		}
		if(max_prev_col[jj].x<score.z) {
			max_prev_col[jj].x=score.z;
			cond1=2;
		}

		max_prev_col[jj].y=max_prev_row_i_ii.y=score.x+gap_open;
		if(max_prev_row_i_ii.y<score.y) {
			max_prev_row_i_ii.y=score.y;
			cond2=1;
		}
		auto gap_open_y=score.y+gap_open;
		if(max_prev_col[jj].y<gap_open_y) {
			max_prev_col[jj].y=gap_open_y;
			cond3=1;
		}
		auto gap_open_z=score.z+gap_open;
		if(max_prev_row_i_ii.y<gap_open_z) {
			max_prev_row_i_ii.y=gap_open_z;
			cond2=2;
		}
		if(max_prev_col[jj].y<score.z) {
			max_prev_col[jj].y=score.z;
			cond3=2;
		}

		if(debug)
		printf("@%u score: %u %u (%u %u %d): %d %d %d\n", thr_id, i+ii, j+jj, a_cache[ii].value, b_cache[jj].value, bp_score, score.x);
		if(/*i+ii<a_size && */score.x>_top.score) {
			_top={i+ii, j+jj, score.x};
		}
		auto& tb_tmp=tb_batch[jj/8+(tile_size.y/8)*ii];
		tb_tmp=(tb_tmp<<4)|encode_conds(cond1, cond2, cond3);
	}
	max_prev_row[i+ii]=max_prev_row_i_ii;
}
			auto ptr=&mat2[j*max_query_len/8+i*(tile_size.y/8)];
			for(unsigned int jj=0; jj<tile_size.y/8*tile_size.x; ++jj) {
				ptr[jj]=reinterpret_cast<unsigned int*>(tb_batch)[jj];
			}
		}
	}

	auto tb_offset=cigar_s[work_id];
	unsigned int tb_size=cigar_s[work_id+1]-tb_offset;
	//assert(tb_size>=a_size/4+1);
	//assert(tb_size>=48+1);
	auto tb_ops=&cigar_buf[tb_offset];

	unsigned int tb_i=0;
	unsigned int i=_top.ii;
	unsigned int j=_top.jj;

#if 1
	enum {
		Mat,
		Ins,
		Del,
		End,
	} st=Mat;
	int s=1;
	int sss=_top.score;
	unsigned int tmp_cnt{0};
	unsigned char tmp_op{0b100000};
	unsigned int tmp_cache{0};
	auto add_op=[&tmp_cnt,&tmp_op,&tb_i,&tmp_cache,tb_ops,tb_size](unsigned char op) {
		if(op==tmp_op) {
			++tmp_cnt;
		} else {
			constexpr unsigned int per_slot=(0xffu>>2)+1;
			while(tmp_cnt>0) {
				unsigned int n=min(tmp_cnt, per_slot);
				if(tb_i/4+1<tb_size) {
					auto rem=tb_i%4;
					tmp_cache|=(((n-1)<<2)|tmp_op)<<(rem*8);
					if(rem==3) {
						tb_ops[tb_i/4+1]=tmp_cache;
						tmp_cache=0;
					}
					++tb_i;
				}
				tmp_cnt-=n;
			}
			tmp_cnt=1;
			tmp_op=op;
		}
	};
	auto elem2=[mat2,max_query_len](unsigned int i, unsigned int j) ->unsigned int {
		auto jj=j%tile_size.y;
		j-=jj;
		return (mat2[(j*max_query_len+i*tile_size.y+jj)/8]>>((7-jj%8)*4))&0xf;
	};
	do {
		if(debug)
		printf("@%u %u %u %u %u %d\n", thr_id, i, j, eccl::get_code<Type>(query_buf, i).value, eccl::get_code<Type>(target_buf, j).value, sss);
		unsigned int cond;
		switch(st) {
			case Mat:
				{
					auto a=eccl::get_code<Type>(query_buf,
							a_offset+i);
					auto b=eccl::get_code<Type>(target_buf,
							b_offset+j);
					s=subst_mat[eccl::pair<Type>{a, b}.value];
				}
				if(s>0) {
					add_op(0b11);
					//assert(ss>0);
				} else {
					add_op(0b00);
					//assert(ss<=0);
				}
				sss-=s;
				if(sss==0) {
					st=End;
					break;
				}
				assert(i>0);
				assert(j>0);
				cond=decode_conds(elem2(i-1, j-1))/9;
				switch(cond) {
					case 0:
						break;
					case 1:
						st=Del;
						break;
					case 2:
						st=Ins;
						break;
					default:
						assert(0);
				}
				--i;
				--j;
				break;
			case Del:
				add_op(0b10);
				sss-=gap_ext;
				cond=j>0?decode_conds(elem2(i, j-1))/3%3:0;
				switch(cond) {
					case 0:
						sss-=gap_open;
						st=Mat;
						break;
					case 1:
						break;
					case 2:
						sss-=gap_open;
						st=Ins;
						break;
					default:
						assert(0);
				}
				--j;
				break;
			case Ins:
				add_op(0b01);
				sss-=gap_ext;
				cond=i>0?decode_conds(elem2(i-1, j))%3:0;
				switch(cond) {
					case 0:
						sss-=gap_open;
						st=Mat;
						break;
					case 1:
						sss-=gap_open;
						st=Del;
						break;
					case 2:
						break;
					default:
						assert(0);
				}
				--i;
				break;
			case End:
				assert(0);
		}
	} while(st!=End);
	add_op(0b100000);
	if((tb_i%4)>0)
		tb_ops[tb_i/4+1]=tmp_cache;
	assert((tb_i+3)/4<tb_size);
#endif
	_top.i0=i;
	_top.j0=j;
	tb_ops[0]=tb_i;

	auto& output=results[work_id];
	output.score=_top.score;
	output.query_se={_top.i0, _top.ii};
	output.target_se={_top.j0, _top.jj};
}

template<eccl::seq_type Type>
void eccl::sw_gpu(const sw_gpu_opts& opts, sw_gpu_args&& args) {

	align_kernel<Type><<<opts.num_blocks, opts.num_threads, 0, args.stream>>>(
			args.work_size, args.max_query_len,
			opts.band_size,
			opts.subst_mat,
			args.target_buf, args.target_se,
			args.query_buf, args.query_se,
			args.tmp_scores,
			args.tmp_traceback_buf, args.tmp_traceback_s,
			args.results,
			args.cigar_buf, args.cigar_s,
			0);

}

template void eccl::sw_gpu<eccl::seq_type::xna>(const sw_gpu_opts& opts, sw_gpu_args&& args);
template void eccl::sw_gpu<eccl::seq_type::prot>(const sw_gpu_opts& opts, sw_gpu_args&& args);

