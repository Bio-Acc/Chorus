#include "qit.h"

QIT::QIT()
{
    // init
    index_size = (uint8_t *)calloc(qit_length, sizeof(uint8_t));
    q_num = (int *)calloc(qit_length * qit_width, sizeof(int));
    q_idx = (int *)calloc(qit_length * qit_width, sizeof(int));
}

QIT::~QIT()
{
    free(index_size);
    free(q_num);
    free(q_idx);
}

int QIT::get_index_size(string mers)
{
    return index_size[get_mers_num(mers)];
}

int QIT::find_empty_dummy(int loc)
{
    int pw = 1;
    int pwn = 32;
    int dummy = loc + 26 * pw - loc % pwn;
    while (dummy < qit_length)
    {
        int area_size = (32 - 26) * pw;
        for (int i = 0; i < area_size; i++)
            if (index_size[dummy + i] == 0)
            {
                return dummy + i;
            }
        pw *= 32;
        pwn *= 32;
        dummy = loc + 26 * pw - loc % pwn;
    }
    // int dummy = loc;

    // dummy = dummy + 26 - dummy % 32;
    // while (dummy % 32 < 30)
    // {
    //     if (this->index_size[dummy] == 0)
    //         break;
    //     dummy++;
    // }
    // if (dummy % 32 != 0)
    // {
    //     return dummy;
    // }

    // dummy = dummy + 26 * 32 - dummy % 1024;
    // while (dummy % 1024 < 960)
    // {
    //     if (this->index_size[dummy] == 0)
    //         break;
    //     dummy++;
    // }
    // if (dummy % 1024 != 0)
    // {
    //     return dummy;
    // }
    printf("out of qit size\n");
    exit(-1);
}

void QIT::add(string mers, int num, int idx)
{
    bool same_f = FILTER_ALL_SAME_AA_SEED;
    for (int i = 0; i < seed_length; i++)
    {
        if (mers[i] < 0)
            return;
        if (i > 0 && mers[i] != mers[i - 1])
            same_f = false;
    }
    if (same_f)
        return;

    int loc = get_mers_num(mers);
    int b_loc = loc;

    int offset = loc * qit_width;
    
    index_size[loc]++;
    if (index_size[loc] <= qit_width)
    {
        q_num[offset + index_size[loc] -1 ] = num;
        q_idx[offset + index_size[loc] -1 ] = idx;
        return;
    }
    

    while (q_num[offset + qit_width - 1] == -1)
    {
        loc = loc + q_idx[offset + qit_width - 1];
        offset = loc * qit_width;
    }
    if (q_num[offset + qit_width - 1] != 0)
    {
        int last_q_num = q_num[offset + qit_width - 1];
        int last_q_idx = q_idx[offset + qit_width - 1];
        q_num[offset + qit_width - 1] = -1;
        int dummy_loc = find_empty_dummy(b_loc);
        q_idx[offset + qit_width - 1] = dummy_loc - loc;

        loc = dummy_loc;
        offset = loc * qit_width;
        q_num[offset] = last_q_num;
        q_idx[offset] = last_q_idx;
        index_size[loc]++;
    }

    int size = index_size[loc];
    q_num[offset + size] = num;
    q_idx[offset + size] = idx;
    index_size[loc]++;
    
}

// void QIT::visualize(int begin, int end)
// {
//     for (int i = begin; i < end; i++)
//     {
//         cout << ITOC[i / 1024] << ITOC[(i % 1024) / 32] << ITOC[i % 32] << " " << index_size[i] << " | ";
//         for (int j = 0; j < WIDTH; j++)
//         {
//             cout << q_num[i * WIDTH + j] << " ";
//         }
//         cout << "| ";
//         for (int j = 0; j < WIDTH; j++)
//         {
//             cout << q_idx[i * WIDTH + j] << " ";
//         }
//         cout << endl;
//     }
// }