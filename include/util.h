#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cmath>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <assert.h>
#include <dirent.h>
#include <limits>
#include <sys/time.h>
#include <string>
#include <bitset>
#include <cstring>

#include <algorithm>
#include <stdio.h>
#include <ctime>
#include "params.h"
#include "qit.h"
#include <thread>
#include <mutex>
#include <atomic>
#include "ThreadPool.h"
#include <chrono>
#include <regex>
#include "omp.h"
// #include "boost/asio.hpp"

extern ThreadPool* pool;
// extern StripedSmithWaterman::Aligner aligner;
// extern StripedSmithWaterman::Filter filter;
// typedef boost::asio::thread_pool ThreadPool;

#define CUDA_CALL(F)                                                          \
    if ((F) != cudaSuccess)                                                   \
    {                                                                         \
        printf("Error %s at %s:%d\n", cudaGetErrorString(cudaGetLastError()), \
               __FILE__, __LINE__);                                           \
        exit(-1);                                                             \
    }

using namespace std;

struct QueryGroup
{
    uint32_t group_id;
    uint32_t n_query = 0;
    uint32_t id[MAX_QUERY_PER_GROUP];
    uint32_t offset[MAX_QUERY_PER_GROUP];
    uint32_t length[MAX_QUERY_PER_GROUP];
    uint32_t hashtable_capacity[MAX_QUERY_PER_GROUP];
    uint32_t hashtable_offset[MAX_QUERY_PER_GROUP];
    uint32_t min_diag_hit[MAX_QUERY_PER_GROUP];
    QIT qit;
};

struct AlignInfo
{
    int group_id;
    int stream_id;
    int idx;
};

struct SWTasks
{
    const char *q;
    const char *c;
    char* c_all[NUM_STREAM];
    size_t c_offs[NUM_STREAM];
    size_t num_task = 0;
    size_t c_len = 0;
    size_t c_offset = 0;
    vector<uint32_t> q_idxs;
    vector<uint32_t> q_lens; 
    // vector<size_t> s_idxs;
    // vector<size_t> s_lens;
    vector<size_t> diags;
    // vector<size_t> q_len4_sum;

    // vector<size_t> s_lim_begin;
    // vector<size_t> s_lim_end;
    vector<size_t> q_len4_offs;
    vector<size_t> s_len4_offs;
    // size_t total_q_len4 = 0;
    // size_t total_s_len4 = 0;

    vector<AlignInfo> info;
};

struct SWResult
{
    bool report = true;
    int num_q;
    // size_t num_s;
    vector<int> q_res;
    vector<size_t> s_res;
    int score;
    int bitscore;
    int align_length;
    double p_identity;
    int n_identity;
    int mismatch;
    int positive;
    int gap_open;
    int gaps;
    string q;
    string s;
    string s_ori;
    string match;
    size_t begin_q;
    size_t end_q;
    size_t begin_s;
    size_t end_s;
    double e_value;

    size_t sn_offset;
    size_t sn_len;
    uint32_t s_len;
    string s_name;

    SWResult(int q)
    {
        num_q = q;
    }
    SWResult()
    {
        
    }
};

// class DisplayResult
// {
// public:
//     // int num_q;
//     // string q;
//     // string s;
//     // string s_ori;
//     // size_t begin_q;
//     // size_t end_q;
//     // size_t begin_s;
//     // size_t end_s;
//     // int align_length;
//     // int score;
//     // int bitscore;
//     // double p_identity;
//     // int n_identity;
//     // int mismatch;
//     // int positive;
//     // int gap_open;
//     // int gaps;

//     // uint32_t db_id;
//     // uint32_t num_stream;
//     // double e_value;
//     size_t sn_offset;
//     size_t sn_len;
//     uint32_t s_len;
    
//     string s_name;

//     DisplayResult()
//     {
        
//     }

//     DisplayResult(SWResult sw)
//     {
//         num_q = sw.num_q;
//         score = sw.score;
//         // num_s = sw.num_s;
//         bitscore = sw.bitscore;
//         align_length = sw.align_length;
//         p_identity = sw.p_identity;
//         n_identity = sw.n_identity;
//         mismatch = sw.mismatch;
//         positive = sw.positive;
//         gap_open = sw.gap_open;
//         gaps = sw.gaps;
//         q = sw.q;
//         s = sw.s;
//         s_ori = sw.s_ori;
//     }
// };

double timeuse(struct timeval start_time, struct timeval end_time);
int get_pos(size_t idx, const size_t *offsets, size_t n_subj);
uint32_t load_fasta(const char *file, char *&str, vector<uint32_t> &offsets, vector<string> &names);
size_t load_offsets(string db, size_t *&s_offsets, size_t *&n_offsets);
int check_db(const char *db, size_t &max_size, size_t& total_size);
// string get_name(string db_name, size_t offset, size_t len);
void load_seq(string db, int num, char *&str, size_t &len);
void proceed_result(vector<SWResult> *res_d, vector<SWResult> &res_t, const char *query, const char *subj, QueryGroup &q_group, const char *s_name, const size_t* s_offsets, const size_t* sn_offsets, const size_t s_num ,size_t total_db_size);

char get_char(const char *s, size_t offset);

void get_arg(const char* name, int& v, int default_v);
void get_arg(const char* name, double& v ,double default_v);
void get_arg(const char* name, string& v);
void get_arg(const char *name, vector<string> &v);

void load_must_include();
bool check_include(const char *s, size_t begin_p, size_t end_p);
bool check_include(string str);