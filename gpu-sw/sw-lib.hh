#ifndef GLF_SW_LIB_HH
#define GLF_SW_LIB_HH

#include <vector>
#include <array>
#include <optional>

//refactor!!!
#include <cstdint>
#include <vector>
#include <string>
#include <string_view>

#include "core.hh"

#define MAPSTR "MXDIE"
#ifndef CIGAR_SHIFT
#define CIGAR_SHIFT 2u
#endif

#define MAX_NUM_CIGAR 16

typedef struct {
    uint16_t sw_score;
    uint16_t mismatch;
    int32_t ref_begin;
    int32_t ref_end;
    int32_t query_begin;
    int32_t query_end;
    uint32_t cigar[MAX_NUM_CIGAR];
    int32_t cigar_len;
} sw_align;

extern const uint8_t encoded_ops[];

static inline uint32_t to_cigar_int(uint32_t length, char op) {
    return (length << CIGAR_SHIFT) | (encoded_ops[(int)op]);
}

static inline char cigar_int_to_op(uint32_t cigar_int) {
    return MAPSTR[cigar_int & 0x3U];
}

static inline uint32_t cigar_int_to_len(uint32_t cigar_int) {
    return cigar_int >> CIGAR_SHIFT;
}


// static void cigar_to_string(const sw_align& alignment, std::string_view& query, std::string_view& ref, std::vector<int>& q_res, std::vector<size_t>& r_res);
static void cigar_to_string(const sw_align& alignment, std::string_view& query, std::string_view& ref, std::vector<int>& q_res, std::vector<size_t>& r_res) {
  q_res.clear();
  r_res.clear();
  int qp = alignment.query_begin;
  int rp = alignment.ref_begin;
  for (int c=0;c<alignment.cigar_len;c++)
  {
    int cigar_int = cigar_int_to_len(alignment.cigar[c]);
    char cigar_op = cigar_int_to_op(alignment.cigar[c]);
    // cout << cigar_int << cigar_op<<endl;
    if (cigar_op=='M'||cigar_op=='X'){
      // q_res.append(query+qp, query+qp+cigar_int);
      for (int i=0;i<cigar_int;i++)
      {
        q_res.push_back(qp+i);
        r_res.push_back(rp+i);
      }
      // r_res.append(ref+rp, ref+rp+cigar_int);
      qp+=cigar_int;
      rp+=cigar_int;
    }
    else if (cigar_op=='I')
    {
      // q_res.append(query+qp, query+qp+cigar_int);
      for (int i=0;i<cigar_int;i++)
      {
        q_res.push_back(qp+i);
        r_res.push_back(-1);
      }
      qp+=cigar_int;
    }
    else if (cigar_op=='D')
    {
      // r_res.append(ref+rp, ref+rp+cigar_int);
      for (int i=0;i<cigar_int;i++)
      {
        r_res.push_back(rp+i);
        q_res.push_back(-1);
      }
      rp+=cigar_int;
    }
    else if (cigar_op=='S')
    {
      // qp+=cigar_int;
    }
    else{
      // std::cout <<"error cigar:"<<cigar_op<<std::endl;
      exit(-1);
    }
  }
}
//!!!

namespace glf {

struct sw_opts {
	struct {
		std::optional<int> match;
		std::optional<int> mismatch;
		std::optional<int> gap_open;
		std::optional<int> gap_extend;
	} scores;
};

struct sw_batch_opts {
	//minscore
};

struct sw_stats {
	std::size_t dev_mem_usage;
	std::size_t host_mem_usage;
};

template<eccl::seq_type T>
struct sw_handle;

template<eccl::seq_type T>
sw_handle<T>* sw_create(const sw_opts& opts);

template<eccl::seq_type T>
sw_handle<T>* sw_create();

template<eccl::seq_type T>
void sw_destroy(sw_handle<T>* hdl);

template<eccl::seq_type T>
std::vector<sw_align> sw_batch(sw_handle<T>* hdl,
		const std::vector<std::string_view>& query,
		const std::vector<std::string_view>& target,
		const sw_batch_opts& opts);

template<eccl::seq_type T>
sw_stats sw_inspect(sw_handle<T>* hdl);

}

#endif
