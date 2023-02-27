#include "query_group.h"

__inline__ size_t min_full_bit(size_t x)
{
    x = x - 1;         // 0010 1100 0000 0000 0000 0000 0000 0000 0000 0001
    x = x | (x >> 1);  // 0011 1110 0000 0000 0000 0000 0000 0000 0000 0000
    x = x | (x >> 2);  // 0011 1111 1000 0000 0000 0000 0000 0000 0000 0000
    x = x | (x >> 4);  // 0011 1111 1111 1000 0000 0000 0000 0000 0000 0000
    x = x | (x >> 8);  // 0011 1111 1111 1111 1111 1000 0000 0000 0000 0000
    x = x | (x >> 16); // 0011 1111 1111 1111 1111 1111 1111 1111 1111 1111
    return (x + 1);    // 0100 0000 0000 0000 0000 0000 0000 0000 0000 0000
}

vector<QueryGroup> init_query_group(int n_query, size_t db_size, vector<uint32_t> q_lengths, vector<uint32_t> q_offsets, char *query, size_t &max_hashtable_capacity, uint32_t &max_n_query)
{

    assert(n_query == q_lengths.size());
    size_t free_byte, total_byte;
    CUDA_CALL(cudaMemGetInfo(&free_byte, &total_byte));
    size_t empty_byte = 536870912;
    empty_byte += qit_length * qit_width * sizeof(int) * MAX_GROUPS_PER_ROUND * 2; 
    empty_byte += qit_length * sizeof(uint8_t) * MAX_GROUPS_PER_ROUND;
    
    empty_byte += MAX_FILTER_TASK * sizeof(Task) * NUM_STREAM;
    empty_byte += sizeof(uint32_t) * NUM_STREAM;

    //  8 * MAX_GROUPS_PER_ROUND * QIT_LENGTH * QIT_WIDTH * sizeof(int)
    free_byte = (free_byte - db_size - empty_byte) / (NUM_STREAM) / sizeof(KeyValue);
    // free_byte = free_byte > numeric_limits<uint32_t>::max() ? numeric_limits<uint32_t>::max() : free_byte;
    // free_byte = 5000000;
    cout << "Free byte = " << (double)free_byte * sizeof(KeyValue) * NUM_STREAM / 1024 / 1024 / 1024 << " GB" << endl;
    uint32_t cap[n_query];
    int id[n_query];

    int hashtable_size_ratio;
    get_arg("hash-size", hashtable_size_ratio, D_HASH_RATIO);
    for (int i = 0; i < n_query; i++)
    {
        cap[i] = (min_full_bit(((db_size * (q_lengths[i]<128?128:q_lengths[i]) / (int)pow(10, seed_length))) / NUM_STREAM) << hashtable_size_ratio )>> 2;
        // cap[i] = (min_full_bit(q_lengths[i]) * (min_full_bit(db_size) >> HASHTABLE_SIZE_RATIO) / NUM_STREAM);
        id[i] = i;
    }

    // greedy algorithm: for each time find largest fitable query, if no fitable query, open a new group

    sort(id, id + n_query, [&](const int &id1, const int &id2)
         { return (cap[id1] == cap[id2]) ? (q_lengths[id1] < q_lengths[id2]) : (cap[id1] < cap[id2]); });

    if (cap[id[n_query - 1]] > free_byte)
    {
        printf("Error: Not enough GPU memory. Please reduce the size of each block of database.\n");
        printf("Longest Requirement %f GB, length = %u, free mem %F GB\n", (double)cap[id[n_query - 1]] / 1024 / 1024 / 1024 * sizeof(KeyValue) * NUM_STREAM, q_lengths[id[n_query - 1]], (double)free_byte / 1024 / 1024 / 1024 * sizeof(KeyValue) * NUM_STREAM);
        exit(1);
    }

    vector<QueryGroup> groups;
    groups.resize(MAX_NUM_GROUPS);
    int n_groups = 0;
    int allocated = 0;
    max_hashtable_capacity = 0;
    max_n_query = 0;
    uint32_t max_len_query = 0;
    while (allocated < n_query)
    {
        n_groups++;
        if (n_groups >= MAX_NUM_GROUPS)
        {
            printf("Max num query groups exceeded!\n");
            exit(-1);
        }
        QueryGroup &g = groups[n_groups - 1];
        g.n_query = 0;
        uint32_t total_hashtable_capacity = 0;
        for (int i = 0; i < n_query; i++)
        {

            if (id[i] == -1)
                continue;
            if (total_hashtable_capacity + cap[id[i]] <= free_byte)
            {
                g.id[g.n_query] = id[i];
                g.hashtable_capacity[g.n_query] = cap[id[i]];
                g.hashtable_offset[g.n_query] = total_hashtable_capacity;
                g.n_query++;

                total_hashtable_capacity += cap[id[i]];
                id[i] = -1;
                allocated++;
            }
            if (g.n_query > max_n_query)
            {
                max_n_query = g.n_query;
            }
            if (g.n_query >= MAX_QUERY_PER_GROUP)
            {
                printf("Max query per group exceeded!\n");
                // exit(-1);
                break;
            }
        }
        if (total_hashtable_capacity > max_hashtable_capacity)
        {
            max_hashtable_capacity = total_hashtable_capacity;
        }
    }

    cout << "Hashtable Capacity = " << (double)max_hashtable_capacity * sizeof(KeyValue) * NUM_STREAM / 1024 / 1024 / 1024 << " GB" << endl;
    cout << "Query index table size = " << (double)(qit_length * qit_width * sizeof(int) * 2) / (1073741824) << " GB" << endl;
    int sens_level;
    get_arg("filter-level", sens_level, D_FILTER_LEVEL);
    sens_level = sens_level + 5 - seed_length;
    for (int i = 0; i < n_groups; i++)
    {
        // printf("Query group %d:[ ", i);
        QueryGroup &g = groups[i];
        g.group_id = i;
        for (int j = 0; j < g.n_query; j++)
        {
            uint32_t id = g.id[j];
            g.length[j] = q_lengths[id];
            g.offset[j] = q_offsets[id];
            g.min_diag_hit[j] = log10(g.length[j]) + sens_level;
            // g.min_diag_hit[j] = MIN_DIAG_HIT;
            // printf("%u(%u,%u,%u) ", g.id[j], g.length[j], g.hashtable_capacity[j], g.min_diag_hit[j]);
            if (g.length[j]>max_len_query)
            {
                max_len_query = g.length[j];
            }
        }
        // printf("]\n");
    }
#pragma omp parallel for
    for (int i = 0; i < n_groups; i++)
    {
        cout << "Building Query Index Table G" << i << ": n_query="<< groups[i].n_query<< "..";
        build_qit(query, groups[i].n_query, groups[i].length, groups[i].offset, groups[i].qit);
    }
    groups.resize(n_groups);
    cout << "Done building query group. Max query length = "<< max_len_query << endl;
    return groups;
}

void build_qit(const char *query, uint32_t n_query, uint32_t *lengths, uint32_t *offsets, QIT &qit)
{
    for (int i = 0; i < n_query; i++)
    {
        int len = lengths[i];
        string str(&query[offsets[i]], &query[offsets[i]] + lengths[i]);
        for (int j = 0; j <= len - seed_length; j++)
        {
            string mer(str.begin() + j, str.begin() + j + seed_length);
            qit.add(mer, i, j);
        }
    }
    for (int i = 0; i < 26; i++)
    {
        qit.index_size[(i << 10) + (i << 5) + i] = 0;
    }

    uint8_t max_mers = 0;
    int nonzero = 0;
    for (int i = 0; i < qit_length; i++)
    {
        if (qit.index_size[i] > 0)
            nonzero++;
        if (qit.index_size[i] > max_mers)
            max_mers = qit.index_size[i];
        // if (qit.index_size[i] > 50)
        //     cout << (char)((i >> 20) + 65) << (char)((i >> 15) % 32 + 65) << (char)((i >> 10) % 32 + 65) << (char)((i >> 5) % 32 + 65) << (char)((i) % 32 + 65) << ":" << (int)qit.index_size[i] << endl;
    }
    // for (int i = 0; i < 26; i++)
    //     for (int j = 0; j < 26; j++)
    //         for (int k = 0; k < 26; k++)
    //         {
    //             int v = qit.index_size[(i << 10) + (j << 5) + k];
    //             if (v > 0)
    //             {
    //                 // cout << (char)(k+65)<< (char)(j+65)<< (char)(i+65)<<v<<endl;
    //                 nonzero++;
    //             }

    //             if (v > max_mers)
    //                 max_mers = v;
    //         }
    cout << "Done. Max element = " << (int)max_mers << ", Nonzero = " << nonzero << " / " << qit_length << endl;
}