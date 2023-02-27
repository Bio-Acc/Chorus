#include "blastp.h"
#include <assert.h>

#define PACK_KEY(k) ((k & ~0x7) | 0x3)

ThreadPool *pool;
mutex mu2;

// vector<SWResult> res_s[MAX_GROUPS_PER_ROUND][NUM_STREAM];

__constant__ uint32_t kHashTableCapacity_dev[MAX_GROUPS_PER_ROUND][MAX_QUERY_PER_GROUP];
__constant__ uint32_t kHashTableOffset_dev[MAX_GROUPS_PER_ROUND][MAX_QUERY_PER_GROUP];

__constant__ int SEED_LENGTH;
__constant__ int QIT_WIDTH;
__constant__ uint32_t MASK;

// 32 bit Murmur3 hash
inline __device__ uint32_t my_hash(uint32_t k, uint32_t kHashTableCapacity)
{
    k ^= k >> 16;
    k *= 0x85ebca6b;
    k ^= k >> 13;
    k *= 0xc2b2ae35;
    k ^= k >> 16;
    return k & (kHashTableCapacity - 1);
}

__device__ int insert_ot(KeyValue *hashtable, uint32_t kHashTableCapacity, uint32_t key, uint32_t value)
{
    key = PACK_KEY(key);
    uint32_t slot = my_hash(key, kHashTableCapacity);
    uint32_t b_slot = slot;
    while (true)
    {
        uint32_t prev = atomicCAS(&hashtable[slot].key, kEmpty, key);
        if (prev == kEmpty || prev == key)
        {
            // hashtable[slot].value = value;
            atomicAdd(&hashtable[slot].value, value);
            return 0;
        }
        slot = (slot + 1) & (kHashTableCapacity - 1);
        if (slot == b_slot)
        {
            return -1;
        }
    }
}

__global__ void seeding_kernel(KeyValue *ht, uint32_t *subj, size_t s_length_block, size_t s_length_total, const uint32_t *q_lengths, const int *q_num, const int *q_idx, int n_query, uint8_t *index_size_dev, uint32_t group_id)
{
    size_t s_begin = ((blockIdx.x * blockDim.x + threadIdx.x) * s_length_block) * 32;

    size_t s_len = s_length_block * 32;
    if (s_begin + s_len >= s_length_total - SEED_LENGTH)
        s_len = s_length_total - SEED_LENGTH - s_begin;
    if (s_len <= 0)
        return;

    size_t s_end = s_begin + s_len;

    for (size_t i = s_begin; i < s_end; i++)
    {
        size_t n_bit = i * 5;
        size_t pos = (n_bit >> 5);
        uint32_t mod = n_bit & 31;
        // assert(pos % 4 == 0);
        uint32_t qit_idx = (subj[pos] >> mod) & MASK;
        if (mod > (31 - (5 * SEED_LENGTH)))
        {
            qit_idx |= (subj[pos + 1] << (32 - mod)) & MASK;
        }

        int hit_size = index_size_dev[qit_idx];

        if (hit_size <= 0)
            continue;

        int qit_p = 0;
        for (int j = 0; j < hit_size; j++)
        {
            int pos = qit_idx * QIT_WIDTH + qit_p;
            int q_num_now = q_num[pos];
            int q_idx_now = q_idx[pos];
            if (q_num_now == -1)
            {
                qit_idx += q_idx_now;
                qit_p = 0;
                pos = qit_idx * QIT_WIDTH;
                q_num_now = q_num[pos];
                q_idx_now = q_idx[pos];
            }

            // printf("%d %d\n",q_num[qit_idx*qit_width+qit_p],q_idx[qit_idx*qit_width+qit_p]);
            unsigned int diag = q_lengths[q_num_now] + i - q_idx_now;
            // KeyValue *pHashTable_addr = ot + q_num_now * kHashTableCapacity_dev[q_num_now];
            KeyValue *pHashTable_addr = ht + kHashTableOffset_dev[group_id][q_num_now];
            int err = insert_ot(pHashTable_addr, kHashTableCapacity_dev[group_id][q_num_now], diag, 1);
            // assert(err != -1);
            if (err == -1)
            {
                printf("Voting Hash Table Full! G%uQ%uK%u\n", group_id, q_num_now, kHashTableCapacity_dev[group_id][q_num_now]);
            }
            qit_p++;
        }
    }
}

__global__ void filter_kernel(KeyValue *ht, Task *tasks, uint32_t *num_task, uint32_t *threshold, uint32_t group_id)
{
    uint32_t q_id = blockIdx.x;
    KeyValue *h_begin = ht + kHashTableOffset_dev[group_id][q_id];

    size_t each_length = (kHashTableCapacity_dev[group_id][q_id] - 1) / blockDim.x + 1;
    h_begin += each_length * threadIdx.x;
    KeyValue *h_end = h_begin + each_length;

    KeyValue *total_end = ht + kHashTableOffset_dev[group_id][q_id] + kHashTableCapacity_dev[group_id][q_id];
    h_end = h_end > total_end ? total_end : h_end;

    Task *task_begin = tasks;

    for (KeyValue *kv = h_begin; kv < h_end; kv++)
    {
        if (kv->key != kEmpty && kv->value != kEmpty && kv->value >= threshold[q_id])
        {
            uint32_t idx = atomicAdd(num_task, 1);
            if (idx >= MAX_FILTER_TASK)
            {
                printf("Filter Task Vector Full! G%uQ%uT%u\n", group_id, q_id, idx);
                return;
            }
            task_begin[idx].key = kv->key;
            task_begin[idx].value = kv->value;
            task_begin[idx].q_id = q_id;
        }
    }

    // size_t total_length = kHashTableOffset_dev[group_id][n_query-1] + kHashTableCapacity_dev[group_id][n_query-1];
    // size_t each_length = (total_length-1)/b + 1;
}

#ifdef USE_GPU_SW
void handle_results(cudaEvent_t &stream, Task *task_host, uint32_t *num_task, QueryGroup &q_group, size_t s_length, int stream_id, vector<SWResult> &res, SWTasks &sw_task)
{
    cudaEventSynchronize(stream);
    mu2.lock();
    cout << "=";
    res.clear();
    size_t n_task_pre = sw_task.num_task;
    size_t n_task = sw_task.num_task + *num_task;
    sw_task.c_len += s_length;
    sw_task.q_idxs.resize(n_task);
    sw_task.q_lens.resize(n_task);
    sw_task.q_len4_offs.resize(n_task+1);
    sw_task.s_len4_offs.resize(n_task+1);
    sw_task.diags.resize(n_task);
    sw_task.info.resize(n_task);
    Task *t_begin = task_host;
    sw_task.num_task = n_task;
    res.resize(*num_task);
#pragma omp parallel for
    for (int i = 0; i < *num_task; i++)
    {
        Task &kv = *(t_begin + i);
        sw_task.q_idxs[i + n_task_pre]=q_group.offset[kv.q_id];
        sw_task.q_lens[i + n_task_pre]=q_group.length[kv.q_id];
        sw_task.diags[i + n_task_pre] =  sw_task.c_offset + kv.key;
        sw_task.info[i+ n_task_pre].group_id = q_group.group_id;
        sw_task.info[i+ n_task_pre].stream_id = stream_id;
        sw_task.info[i+ n_task_pre].idx = i;
        res[i].num_q = kv.q_id;
    }

    for (int i = 0; i < *num_task; i++)
    {
        int q_len4 = sw_task.q_lens[i+ n_task_pre];
        q_len4 = q_len4 % 4? q_len4 + (4 - (q_len4 % 4)) : q_len4;
        int s_len4 = sw_task.q_lens[i+ n_task_pre] + (band_width << 1);
        s_len4 = s_len4 %4? s_len4 + (4-(s_len4%4)):s_len4;
        sw_task.q_len4_offs[i+ n_task_pre+1] = sw_task.q_len4_offs[i+ n_task_pre] + q_len4;
        sw_task.s_len4_offs[i+ n_task_pre+1] = sw_task.s_len4_offs[i+ n_task_pre] + s_len4;
    }

    sw_task.c_offset += s_length;
    mu2.unlock();

}
#else
void handle_results(cudaEvent_t &stream, const char *query, const char *subj, Task *task_host, uint32_t *num_task, QueryGroup &q_group, size_t s_length, int stream_id, vector<SWResult> &res, SWTasks &sw_task, ThreadPool *pool, vector<future<int>> &rs)
{
    cudaEventSynchronize(stream);
    cout << "=";
    res.clear();
    res.resize(*num_task);
    sw_task.q = query;
    sw_task.c = subj;
    sw_task.c_len = s_length;
    sw_task.q_idxs.resize(*num_task);
    sw_task.q_lens.resize(*num_task);
    sw_task.diags.resize(*num_task);
    Task *t_begin = task_host;
    sw_task.num_task = *num_task;
#pragma omp parallel for
    for (int i = 0; i < *num_task; i++)
    {
        Task &kv = *(t_begin + i);
        sw_task.q_idxs[i]=q_group.offset[kv.q_id];
        sw_task.q_lens[i]=q_group.length[kv.q_id];
        sw_task.diags[i]=kv.key;
        res[i].num_q = kv.q_id;
    }
    mu2.lock();
    for (int i = 0; i < sw_task.num_task; ++i)
    {
        rs.emplace_back(pool->enqueue([&, i]
                                      {
            smith_waterman_kernel(i,&res[i],&sw_task);
            return i; }));
    }
    mu2.unlock();
}
#endif

void search_db_batch(const char *query, char *subj[], vector<QueryGroup> &q_groups, size_t s_length[], Task *task_host[][NUM_STREAM], uint32_t *task_num_host[][NUM_STREAM], size_t max_hashtable_capacity, uint32_t max_n_query, uint32_t total_len_query, string db_name, uint32_t db_num, vector<SWResult> *res, size_t total_db_size, TimeProfile &time_prof)
{
    struct timeval t_start, t_end, tt_start;

    gettimeofday(&t_start, NULL);

    CUDA_CALL(cudaMemcpyToSymbol(SEED_LENGTH, &seed_length, sizeof(int)));
    CUDA_CALL(cudaMemcpyToSymbol(QIT_WIDTH, &qit_width, sizeof(int)));
    uint32_t mask = (uint32_t)pow(2, 5 * seed_length) - 1;
    CUDA_CALL(cudaMemcpyToSymbol(MASK, &mask, sizeof(uint32_t)));

    size_t sum_s_len = 0;
    for (int i = 0; i < NUM_STREAM; i++)
    {
        sum_s_len += s_length[i];
        assert(s_length[i] % 32 == 0);
    }

    char *subj_dev;
    CUDA_CALL(cudaMalloc((void **)&subj_dev, sum_s_len / 8 * 5));

    int n_groups = q_groups.size();
    if (n_groups > MAX_GROUPS_PER_ROUND)
        n_groups = MAX_GROUPS_PER_ROUND;

    int *q_num_dev[n_groups];
    int *q_idx_dev[n_groups];
    uint32_t *q_lengths_dev[n_groups];
    uint8_t *index_size_dev[n_groups];
    uint32_t *threshold_dev[n_groups];

    uint32_t kHashTableCapacity_host[MAX_GROUPS_PER_ROUND][MAX_QUERY_PER_GROUP];
    uint32_t kHashTableOffset_host[MAX_GROUPS_PER_ROUND][MAX_QUERY_PER_GROUP];

    for (int g = 0; g < n_groups; g++)
    {
        CUDA_CALL(cudaMalloc((void **)&q_num_dev[g], qit_length * qit_width * sizeof(int)));
        CUDA_CALL(cudaMalloc((void **)&q_idx_dev[g], qit_length * qit_width * sizeof(int)));
        CUDA_CALL(cudaMalloc((void **)&q_lengths_dev[g], MAX_QUERY_PER_GROUP * sizeof(uint32_t)));
        CUDA_CALL(cudaMalloc((void **)&index_size_dev[g], qit_length * sizeof(uint8_t)));
        CUDA_CALL(cudaMalloc((void **)&threshold_dev[g], MAX_QUERY_PER_GROUP * sizeof(uint32_t)));
    }
    KeyValue *pHashTable[NUM_STREAM];
    Task *task_dev[NUM_STREAM];
    uint32_t *task_num_dev[NUM_STREAM];
    for (int s = 0; s < NUM_STREAM; s++)
    {
        pHashTable[s] = create_hashtable(max_hashtable_capacity);
        CUDA_CALL(cudaMalloc((void **)&task_dev[s], MAX_FILTER_TASK * sizeof(Task)));
        CUDA_CALL(cudaMemset(task_dev[s], 0, MAX_FILTER_TASK * sizeof(Task)));
        CUDA_CALL(cudaMalloc((void **)&task_num_dev[s], sizeof(uint32_t)));
        CUDA_CALL(cudaMemset(task_num_dev[s], 0, sizeof(uint32_t)));
    }

    char *s_name[NUM_STREAM] = {0};
    size_t *s_offsets[NUM_STREAM] = {0};
    size_t *sn_offsets[NUM_STREAM] = {0};
    size_t s_num[NUM_STREAM] = {0};

    int mingridsize_seeding, mingridsize_filter;
    int threadblocksize_seeding, threadblocksize_filter;
    CUDA_CALL(cudaOccupancyMaxPotentialBlockSize(&mingridsize_seeding, &threadblocksize_seeding, seeding_kernel, 0, 0));
    CUDA_CALL(cudaOccupancyMaxPotentialBlockSize(&mingridsize_filter, &threadblocksize_filter, filter_kernel, 0, 0));

    // cout << "Seeding Block size:" << threadblocksize_seeding <<"," << mingridsize_seeding <<endl;
    // cout << "Filter Block size:" << threadblocksize_filter <<"," << mingridsize_filter <<endl;

    size_t free_byte, total_byte;
    CUDA_CALL(cudaMemGetInfo(&free_byte, &total_byte));
    cout << "GPU mem: " << (double)(total_byte - free_byte) / (1073741824) << " GB / " << (double)total_byte / (1073741824) << " GB" << endl;

#ifndef USE_GPU_SW
    SWTasks sw_tasks[q_groups.size()][NUM_STREAM];
#endif
    SWTasks sw_tasks_total;
    vector<SWResult> res_s[q_groups.size()][NUM_STREAM];

    gettimeofday(&t_end, NULL);
    time_prof.mem_time += timeuse(t_start, t_end);

    int g_begin = 0;
    while (g_begin < q_groups.size())
    {
        sw_tasks_total.c_offset = 0;
        double group_time = 0;
        cout << "Group " << g_begin + 1 << "/" << q_groups.size() << "\t[";
        gettimeofday(&t_start, NULL);
        n_groups = q_groups.size() - g_begin;
        if (n_groups > MAX_GROUPS_PER_ROUND)
            n_groups = MAX_GROUPS_PER_ROUND;
        for (int g = g_begin; g < g_begin + n_groups; g++)
        {
            int g_idx = g - g_begin;
            CUDA_CALL(cudaMemcpy(q_num_dev[g_idx], q_groups[g].qit.q_num, qit_length * qit_width * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(q_idx_dev[g_idx], q_groups[g].qit.q_idx, qit_length * qit_width * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(index_size_dev[g_idx], q_groups[g].qit.index_size, qit_length * sizeof(uint8_t), cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(q_lengths_dev[g_idx], q_groups[g].length, MAX_QUERY_PER_GROUP * sizeof(uint32_t), cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(threshold_dev[g_idx], q_groups[g].min_diag_hit, MAX_QUERY_PER_GROUP * sizeof(uint32_t), cudaMemcpyHostToDevice));
            memcpy(kHashTableCapacity_host[g_idx], q_groups[g].hashtable_capacity, sizeof(uint32_t) * MAX_QUERY_PER_GROUP);
            memcpy(kHashTableOffset_host[g_idx], q_groups[g].hashtable_offset, sizeof(uint32_t) * MAX_QUERY_PER_GROUP);
        }

        CUDA_CALL(cudaMemcpyToSymbol(kHashTableCapacity_dev, kHashTableCapacity_host, sizeof(uint32_t) * MAX_QUERY_PER_GROUP * MAX_GROUPS_PER_ROUND));
        CUDA_CALL(cudaMemcpyToSymbol(kHashTableOffset_dev, kHashTableOffset_host, sizeof(uint32_t) * MAX_QUERY_PER_GROUP * MAX_GROUPS_PER_ROUND));

        cudaStream_t streams[NUM_STREAM];
        thread result_threads[n_groups][NUM_STREAM];

        cudaEvent_t seeding_finished[n_groups][NUM_STREAM];

#ifndef USE_GPU_SW
        vector<future<int>> rs[n_groups][NUM_STREAM];
#endif
        size_t s_begin = 0;

        gettimeofday(&t_end, NULL);
        group_time += timeuse(t_start, t_end);
        time_prof.mem_time += timeuse(t_start, t_end);
        // cout << "Prepare mem and data Time: " << timeuse(t_start, t_end) << endl;

        gettimeofday(&t_start, NULL);

        for (int s = 0; s < NUM_STREAM; s++)
        {
            CUDA_CALL(cudaStreamCreate(&streams[s]));
            // printf("start stream %d\n", s);
            size_t s_length_stream = s_length[s];
            size_t s_length_stream_byte = s_length_stream / 8 * 5;
            size_t s_length_stream_block = s_length_stream / 32 * 5;
            size_t each_length_block = (s_length_stream_block - 1) / (mingridsize_seeding * threadblocksize_seeding) + 1;

            if (g_begin == 0)
            {
                CUDA_CALL(cudaMemcpyAsync(subj_dev + s_begin, subj[s], s_length_stream_byte, cudaMemcpyHostToDevice, streams[s]));
            }
            if (STREAM_SYNC && s > 0)
            {
                CUDA_CALL(cudaStreamSynchronize(streams[s - 1]));
            }
            for (int g = g_begin; g < g_begin + n_groups; g++)
            {
                int g_idx = g - g_begin;
                CUDA_CALL(cudaEventCreate(&seeding_finished[g_idx][s]));
                int n_query = q_groups[g].n_query;
                if (g > 0)
                {
                    CUDA_CALL(cudaMemsetAsync(task_dev[s], 0, MAX_FILTER_TASK * sizeof(Task), streams[s]));
                    CUDA_CALL(cudaMemsetAsync(task_num_dev[s], 0, sizeof(uint32_t), streams[s]));
                    CUDA_CALL(cudaMemsetAsync(pHashTable[s], 0xff, max_hashtable_capacity * sizeof(KeyValue), streams[s]));
                }
                seeding_kernel<<<mingridsize_seeding, threadblocksize_seeding, 0, streams[s]>>>(pHashTable[s], (uint32_t *)(subj_dev + s_begin), each_length_block, s_length_stream, q_lengths_dev[g_idx], q_num_dev[g_idx], q_idx_dev[g_idx], n_query, index_size_dev[g_idx], g_idx);
                filter_kernel<<<n_query, threadblocksize_filter, 0, streams[s]>>>(pHashTable[s], task_dev[s], task_num_dev[s], threshold_dev[g_idx], g_idx);
                // CUDA_CALL(cudaMemcpyAsync(hashtable_host[g_idx][s], pHashTable[s], max_hashtable_capacity * sizeof(KeyValue), cudaMemcpyDeviceToHost, streams[s]));
                CUDA_CALL(cudaMemcpyAsync(task_host[g_idx][s], task_dev[s], MAX_FILTER_TASK * sizeof(Task), cudaMemcpyDeviceToHost, streams[s]));
                CUDA_CALL(cudaMemcpyAsync(task_num_host[g_idx][s], task_num_dev[s], sizeof(uint32_t), cudaMemcpyDeviceToHost, streams[s]));
                CUDA_CALL(cudaEventRecord(seeding_finished[g_idx][s]));
#ifdef USE_GPU_SW
                result_threads[g_idx][s] = thread(handle_results, ref(seeding_finished[g_idx][s]), task_host[g_idx][s], task_num_host[g_idx][s], ref(q_groups[g]), s_length[s], s, ref(res_s[g][s]), ref(sw_tasks_total));
#else
                result_threads[g_idx][s] = thread(handle_results, ref(seeding_finished[g_idx][s]), query, subj[s], task_host[g_idx][s], task_num_host[g_idx][s], ref(q_groups[g]), s_length[s], s, ref(res_s[g][s]), ref(sw_tasks[g][s]), pool, ref(rs[g_idx][s]));
#endif 
            }
            s_begin += s_length_stream_byte;
            cout << "=";
        }

        CUDA_CALL(cudaDeviceSynchronize());

        gettimeofday(&t_end, NULL);
        time_prof.gpu_time += timeuse(t_start, t_end);
        group_time += timeuse(t_start, t_end);
        // cout << "GPU computing Time: " << timeuse(t_start, t_end) << endl;

        if (g_begin == 0)
        {
            gettimeofday(&tt_start, NULL);
            for (int s = 0; s < NUM_STREAM; s++)
            {
                string fname = db_name + "_" + to_string(db_num) + "_" + to_string(s) + ".name";
                int fd = open(fname.data(), O_RDONLY);
                if (fd == -1)
                {
                    std::cerr << "Error opening '" << fname << ". Bailing out." << std::endl;
                    exit(1);
                }
                size_t len = lseek(fd, 0, SEEK_END);
                char *map = (char *)mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);
                close(fd);
                s_name[s] = (char *)malloc(len);
                memcpy(s_name[s], map, len);
                munmap(map, len);

                s_num[s] = load_offsets(db_name + "_" + to_string(db_num) + "_" + to_string(s), s_offsets[s], sn_offsets[s]);
            }
            gettimeofday(&t_end, NULL);
            time_prof.name_time += timeuse(tt_start, t_end);
            group_time += timeuse(tt_start, t_end);
            // cout << "Load seqs name Time: " << timeuse(tt_start, t_end) << endl;
        }

        gettimeofday(&tt_start, NULL);

        for (int s = 0; s < NUM_STREAM; s++)
        {
            CUDA_CALL(cudaStreamDestroy(streams[s]));
        }

        int hsp_count = 0;

        for (int s = 0; s < NUM_STREAM; s++)
        {
            for (int g = g_begin; g < g_begin + n_groups; g++)
            {
                int g_idx = g - g_begin;
                result_threads[g_idx][s].join();
                CUDA_CALL(cudaEventDestroy(seeding_finished[g_idx][s]));
                hsp_count += res_s[g][s].size();
                cout << "=";
#ifndef USE_GPU_SW
                for (auto &r : rs[g_idx][s])
                    r.get();
                proceed_result(res, res_s[g][s], query, subj[s], q_groups[g], s_name[s], s_offsets[s], sn_offsets[s], s_num[s], total_db_size);
                cout << "=";
#endif
            }
        }

        g_begin += MAX_GROUPS_PER_ROUND;

        gettimeofday(&t_end, NULL);
        time_prof.cpu_time += timeuse(tt_start, t_end);
        group_time += timeuse(tt_start, t_end);
        cout << "] " << group_time << "s, " << hsp_count << " HSPs" << endl;
    }

    gettimeofday(&t_start, NULL);

    n_groups = q_groups.size();
    if (n_groups > MAX_GROUPS_PER_ROUND)
        n_groups = MAX_GROUPS_PER_ROUND;

    for (int s = 0; s < NUM_STREAM; s++)
    {
        destroy_hashtable(pHashTable[s]);
        CUDA_CALL(cudaFree(task_dev[s]));
        CUDA_CALL(cudaFree(task_num_dev[s]));
    }

    for (int g = 0; g < n_groups; g++)
    {
        CUDA_CALL(cudaFree(q_num_dev[g]));

        CUDA_CALL(cudaFree(q_idx_dev[g]));

        CUDA_CALL(cudaFree(q_lengths_dev[g]));

        CUDA_CALL(cudaFree(index_size_dev[g]));

        CUDA_CALL(cudaFree(threshold_dev[g]));
    }

    gettimeofday(&t_end, NULL);
    time_prof.mem_time += timeuse(t_start, t_end);

#ifdef USE_GPU_SW
    gettimeofday(&t_start, NULL);
    char* query_dev;
    CUDA_CALL(cudaMalloc((void **)&query_dev, total_len_query));
    CUDA_CALL(cudaMemcpy(query_dev, query, total_len_query, cudaMemcpyHostToDevice));
    sw_tasks_total.q = query;
    for (int s = 0; s < NUM_STREAM; s++)
    {
        sw_tasks_total.c_all[s] = subj[s];
        sw_tasks_total.c_offs[s] = s==0? 0: sw_tasks_total.c_offs[s-1] +s_length[s-1];
    }
    gasal_run(sw_tasks_total, res_s, query_dev, subj_dev, q_groups.size(), band_width);
    cout << "Done.\t[";

    gettimeofday(&t_end, NULL);
    time_prof.gpu_time += timeuse(t_start, t_end);
    gettimeofday(&t_start, NULL);

    CUDA_CALL(cudaFree(query_dev));
    for (int s = 0; s < NUM_STREAM; s++)
    {
        for (int g = 0; g < q_groups.size(); g++)
        {
            proceed_result(res, res_s[g][s], query, subj[s], q_groups[g], s_name[s], s_offsets[s], sn_offsets[s], s_num[s], total_db_size);
        }
        cout << "=";
    }
    cout << "] ";
    gettimeofday(&t_end, NULL);
    cout << timeuse(t_start, t_end) <<"s" << endl;
    time_prof.cpu_time += timeuse(t_start, t_end);
#endif

    gettimeofday(&t_start, NULL);

    CUDA_CALL(cudaFree(subj_dev));

    for (int s = 0; s < NUM_STREAM; s++)
    {
        free(s_name[s]);
        free(sn_offsets[s]);
        free(s_offsets[s]);
    }

    gettimeofday(&t_end, NULL);
    time_prof.mem_time += timeuse(t_start, t_end);
}

void blastp(string argv_query, vector<string> argv_dbs, string argv_out)
{
    vector<uint32_t> q_offsets;
    vector<string> q_names;
    char *query;

    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);
    uint32_t q_length = load_fasta(argv_query.data(), query, q_offsets, q_names);

    q_offsets.push_back(q_length);

    vector<uint32_t> q_lengths;
    for (int i = 0; i < q_offsets.size() - 1; i++)
    {
        q_lengths.push_back(q_offsets[i + 1] - q_offsets[i] - 1);
    }
    int n_query = q_offsets.size() - 1;
    gettimeofday(&t_end, NULL);
    cout << "Load query Time: " << timeuse(t_start, t_end) << endl;
    gettimeofday(&t_start, NULL);

    vector<SWResult> res_d[n_query];

    size_t max_db_size = 0;
    size_t total_db_size = 0;
    vector<int> db_sizes;
    for (int i = 0; i < argv_dbs.size(); i++)
    {
        db_sizes.push_back(check_db(argv_dbs[i].data(), max_db_size, total_db_size));
        if (db_sizes[i] <= 0)
        {
            cout << "DB " << argv_dbs[i] << " not found!" << endl;
            exit(-1);
        }
    }
    cout << "Max db size = " << max_db_size / (1073741824) << " GB" << endl;
    cout << "Total db size = " << (double)total_db_size / (1073741824) << " GB" << endl;
    total_db_size = (total_db_size * 8) / 5;

    size_t max_hashtable_capacity;
    uint32_t max_n_query;
    vector<QueryGroup> q_groups = init_query_group(n_query, max_db_size, q_lengths, q_offsets, query, max_hashtable_capacity, max_n_query);

    // init_hashtable_capacity(n_query, max_db_size, q_lengths);

    // cout << "Each hash table size = ";
    // for (int i = 0; i < n_query; i++)
    // {
    //     cout << (double)kHashTableCapacity_host[i] * sizeof(KeyValue) * NUM_STREAM / (1073741824) << " ";
    // }
    // cout << "GB, total size = " << (double)kHashTableOffset_host[n_query] * sizeof(KeyValue) * NUM_STREAM / (1073741824) << " GB." << endl;

    uint32_t n_groups = q_groups.size();
    if (n_groups > MAX_GROUPS_PER_ROUND)
        n_groups = MAX_GROUPS_PER_ROUND;
    // KeyValue *hashtable_host[n_groups][NUM_STREAM];
    Task *task_host[n_groups][NUM_STREAM];
    uint32_t *task_num_host[n_groups][NUM_STREAM];

    for (int g = 0; g < n_groups; g++)
    {
        for (int s = 0; s < NUM_STREAM; s++)
        {
            // CUDA_CALL(cudaMallocHost(&hashtable_host[g][s], max_hashtable_capacity * sizeof(KeyValue)));
            CUDA_CALL(cudaMallocHost(&task_host[g][s], MAX_FILTER_TASK * sizeof(Task)));
            CUDA_CALL(cudaMallocHost(&task_num_host[g][s], sizeof(uint32_t)));
        }
    }

    pool = new ThreadPool(num_threads);

    gettimeofday(&t_end, NULL);
    cout << "Prepare Time: " << timeuse(t_start, t_end) << endl;

    TimeProfile time_prof;

    struct timeval c_start, c_end;
    gettimeofday(&c_start, NULL);
    for (int d = 0; d < argv_dbs.size(); d++)
    {
        string db_name(argv_dbs[d]);
        for (int i = 0; i < db_sizes[d]; i++)
        {
            struct timeval start, end;
            gettimeofday(&start, NULL);
            cout << "Search DB " << d + 1 << "/" << argv_dbs.size() << ", Part " << i + 1 << "/" << db_sizes[d] << endl;

            char *subj[NUM_STREAM];
            size_t s_size[NUM_STREAM];
            size_t s_len[NUM_STREAM];

            for (int s = 0; s < NUM_STREAM; s++)
            {
                load_seq(db_name + "_" + to_string(i), s, ref(subj[s]), ref(s_size[s]));
                s_len[s] = (s_size[s] * 8) / 5;
            }

            search_db_batch(query, subj, q_groups, s_len, task_host, task_num_host, max_hashtable_capacity, max_n_query, q_length, db_name, i, res_d, total_db_size, time_prof);

            for (int s = 0; s < NUM_STREAM; s++)
            {
                munmap(subj[s], s_size[s]);
            }
            gettimeofday(&end, NULL);
            cout << "Total Batch Time: " << timeuse(start, end) << endl;
        }
    }
    gettimeofday(&c_end, NULL);
    cout << "Finish searching." << endl;
    cout << "GPU Calculation time:\t" << time_prof.gpu_time << endl;
    cout << "CPU Calculation time:\t" << time_prof.cpu_time << endl;
    cout << "Others time:\t" << time_prof.mem_time << endl;
    cout << "Load seqs name Time:\t" << time_prof.name_time << endl;
    cout << "Total Calculation Time:\t" << timeuse(c_start, c_end) << endl;

    gettimeofday(&t_start, NULL);

    for (int g = 0; g < n_groups; g++)
        for (int s = 0; s < NUM_STREAM; s++)
        {
            CUDA_CALL(cudaFreeHost(task_host[g][s]));
            CUDA_CALL(cudaFreeHost(task_num_host[g][s]));
        }

    gettimeofday(&t_end, NULL);
    cout << "Free memory Time:\t" << timeuse(t_start, t_end) << endl;
    gettimeofday(&t_start, NULL);

    for (int i = 0; i < n_query; i++)
    {
        sort_heap(res_d[i].begin(), res_d[i].end(), [&](const SWResult &sw1, const SWResult &sw2)
                  { return (sw1.e_value == sw2.e_value) ? (sw1.score > sw2.score) : (sw1.e_value < sw2.e_value); });
    }

    int outfmt;
    get_arg("outfmt", outfmt, D_OUTFMT);
    switch (outfmt)
    {
    case 0:
        output_result_tabular(argv_out, res_d, query, q_offsets, q_names);
        break;
    case 1:
        output_result_align(argv_out, res_d, query, q_offsets, q_names);
        break;
    case 2:
        output_result_tabular(argv_out, res_d, query, q_offsets, q_names);
        output_result_fa(argv_out + ".fasta", res_d, query, q_offsets, q_names);
        break;
    case 3:
        output_result_cast(argv_out, res_d, query, q_offsets, q_names);
        break;
    case 4:
        output_result_a3m(argv_out, res_d, query, q_offsets, q_names);
        break;
    case 5:
        output_result_reduce(argv_out, res_d, query, q_offsets, q_names);
        break;
    default:
        break;
    }

    free(query);

    gettimeofday(&t_end, NULL);

    cout << "Output Time:\t" << timeuse(t_start, t_end) << endl;

    cout << "Finished." << endl;
}