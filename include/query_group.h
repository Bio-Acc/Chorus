#pragma once

#include "util.h"
#include "qit.h"
#include "hash_table.h"

#if defined(__CUDACC__)
#include "cuda_runtime.h"
#endif

vector<QueryGroup> init_query_group(int n_query, size_t db_size, vector<uint32_t> q_lengths, vector<uint32_t> q_offsets, char *query, size_t &max_hashtable_capacity, uint32_t &max_n_query);

void build_qit(const char *query, uint32_t n_query, uint32_t *lengths, uint32_t *offsets, QIT &qit);
