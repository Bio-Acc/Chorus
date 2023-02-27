#pragma once

#include "util.h"

void output_result_tabular(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names);
void output_result_cast(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names);
void output_result_fa(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names);
void output_result_a3m(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names);
void output_result_reduce(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names);
void output_result_align(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names);

