#pragma once

#include "util.h"

using namespace std;

// void smith_waterman(const char *q, const char *c, const size_t *q_idxs, const size_t *q_lens, const size_t *c_idxs, const size_t *c_lens, size_t num_task, vector<SWResult> &res);

// void banded_smith_waterman(const char *q, const char *c, vector<uint32_t>& q_idxs, vector<uint32_t>& q_lens, vector<size_t>& diags, size_t c_len, size_t num_task, vector<SWResult> &res, ThreadPool* pool, vector<future<int>>& rs);

void smith_waterman_kernel(const int idx, SWResult *res, SWTasks* sw_task);

void gasal_run(SWTasks tasks, vector<SWResult> res[][NUM_STREAM],const char* q_dev, const char* s_dev, int num_g, int span);