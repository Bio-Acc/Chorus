#include "smith.h"
// #include <libunwind.h>
// #include <gperftools/tcmalloc.h>

#define max2(m, n) ((m) > (n) ? (m) : (n))
#define max3(m, n, p) ((m) > (n) ? ((m) > (p) ? (m) : (p)) : ((n) > (p) ? (n) : (p)))
#define MASK5(v) (v & 0b11111)
#define END 0
#define TOP 1
#define LEFT 2
#define DIAG 3

// mutex mu1;

typedef struct
{
    int d;
    int x; // top
    int y; // left
    int m; // left top
    int s;
} record;

inline char get_char(const char *s, size_t offset)
{
    size_t n_bit = offset * 5;
    return MASK5((unsigned)((*((uint16_t *)&(s[n_bit >> 3]))) >> (n_bit & 7)));
}

void smith_waterman_kernel(const int idx, SWResult *res, SWTasks *sw_task)
{
    const char *q = sw_task->q;
    const char *c = sw_task->c;
    size_t c_len = sw_task->c_len;
    size_t q_idx = sw_task->q_idxs[idx];
    size_t n = sw_task->q_lens[idx];
    size_t diag = sw_task->diags[idx]; //  pos of c at end of q

    if (has_must_include)
    {
        int64_t c_begin = (int64_t)diag - band_width - n;
        size_t c_end = diag + band_width;
        c_begin = c_begin >= 0 ? c_begin : 0;
        c_end = c_end <= c_len ? c_end : c_len;
        if (!check_include(c, c_begin, c_end))
        {
            res->report = false;
            return;
        }
    }

    // int64_t c_begin = (int64_t)diag - band_width - n;
    // size_t c_end = diag + band_width;
    // c_begin = c_begin >= 0 ? c_begin : 0;
    // c_end = c_end <= c_len ? c_end : c_len;
    // string sss;
    // printf("%lu,%lu,%lu\n",(size_t)c_begin,c_end,c_len );
    // for (size_t h=c_begin;h<c_end;h++)
    // {
    //     char ch = get_char(c, h)+65;
    //     sss.push_back(ch);
    // }
    // // cout<<sss<<endl;
    

    size_t width = (n + band_width) << 1;
    size_t height = band_width + 2;

    // short *s = (short *)malloc(width * height * sizeof(short));
    // char *p = (char *)malloc(width * height * sizeof(char));
    record *r = (record *)malloc(width * height * sizeof(record));

    if (r == nullptr)
    {
        printf("CPU out of memory!\n");
        exit(-1);
        return;
    }

    // memset(s, 0, width * height * sizeof(short));
    // memset(p, 0, width * height * sizeof(char));
    memset(r, 0, width * height * sizeof(record));

    size_t max_i = 0;
    size_t max_j = 0;

    for (int j = 2; j < width; j++)
    {
        for (int i = 1; i < height - 1; i++)
        {
            int q_pos = i + (j >> 1) - height + 1;
            int64_t c_pos = (int64_t)diag + (j >> 1) + (j & 1) - i - n;
            if (q_pos < 0 || c_pos < 0 || q_pos >= n || c_pos >= c_len)
            {
                // r[j * height + i].s = 0;
                continue;
            }

            // cout<<q_pos<<" "<<c_pos<<endl;
            // s[j * height + i] = 1;

            char chq = q[q_idx + q_pos];
            char chc = get_char(c, c_pos);

            if (chq == END_SIGNAL || chc == END_SIGNAL)
            {
                continue;
            }

            // int score[4] = {0, 0, 0, 0};

            if (j & 1)
            {
                r[j * height + i].x = max3(r[(j - 1) * height + (i - 1)].x + SCORE_GAP_EXT, r[(j - 1) * height + (i - 1)].m + SCORE_GAP, 0);
                r[j * height + i].y = max3(r[(j - 1) * height + i].y + SCORE_GAP_EXT, r[(j - 1) * height + i].m + SCORE_GAP, 0);
            }
            else
            {
                r[j * height + i].x = max3(r[(j - 1) * height + i].x + SCORE_GAP_EXT, r[(j - 1) * height + i].m + SCORE_GAP, 0);
                r[j * height + i].y = max3(r[(j - 1) * height + (i + 1)].y + SCORE_GAP_EXT, r[(j - 1) * height + (i + 1)].m + SCORE_GAP, 0);
            }

            if (chq == ILLEGAL_WORD || chc == ILLEGAL_WORD)
            {
                // illegal word
                r[j * height + i].m = 0;
            }
            else
            {
                r[j * height + i].m = max2(max3(r[(j - 2) * height + i].x, r[(j - 2) * height + i].y, r[(j - 2) * height + i].m) + BLOSUM62[chq * 26 + chc], 0);
            }

            r[j * height + i].s = max3(r[j * height + i].x, r[j * height + i].y, r[j * height + i].m);

            if (r[j * height + i].s != 0)
            {
                if (r[j * height + i].s == r[j * height + i].x)
                    r[j * height + i].d = TOP;
                if (r[j * height + i].s == r[j * height + i].y)
                    r[j * height + i].d = LEFT;
                if (r[j * height + i].s == r[j * height + i].m)
                    r[j * height + i].d = DIAG;
            }

            if (r[j * height + i].s > r[max_j * height + max_i].s)
            {
                max_i = i;
                max_j = j;
            }
        }
    }

    // for (int i = 0; i < height; i++)
    // {
    //     for (int j = 0; j < width; j++)
    //     {
    //         cout << s[j * height + i] << "\t";
    //     }
    //     cout << endl;
    // }

    res->score = r[max_j * height + max_i].s;

    size_t i = max_i;
    size_t j = max_j;
    while (r[j * height + i].s > 0)
    {
        switch (r[j * height + i].d)
        {
        case DIAG:
            res->q_res.push_back(i + (j >> 1) - height + 1 + q_idx);
            res->s_res.push_back(diag + (j >> 1) + (j & 1) - i - n);
            j -= 2;
            break;

        case TOP:
            res->q_res.push_back(i + (j >> 1) - height + 1 + q_idx);
            res->s_res.push_back(-1);
            i -= (j & 1) ? 1 : 0;
            j -= 1;
            break;

        case LEFT:
            res->q_res.push_back(-1);
            res->s_res.push_back(diag + (j >> 1) + (j & 1) - i - n);
            i += (j & 1) ? 0 : 1;
            j -= 1;
            break;

        default:
            printf("err\n");
            break;
        }
    }

    free(r);

    reverse(res->q_res.begin(), res->q_res.end());
    reverse(res->s_res.begin(), res->s_res.end());

    size_t len = res->s_res.size();
    res->align_length = len;
    res->bitscore = (E_lambda * res->score - log(E_k)) / (0.69314718055995);
    char s[len + 1] = {0};
    char q_seq[len + 1] = {0};
    char s_ori[len + 1] = {0};
    char match[len + 1] = {0};
    int s_ori_len = 0;
    res->gap_open = 0;
    res->gaps = 0;
    bool ga = false;
    for (int i = 0; i < len; i++)
    {
        // cout<<c_res[t][i]<<" ";
        if (res->s_res[i] != (size_t)(-1))
        {
            ga = false;
            s[i] = get_char(c, res->s_res[i]) + 65;
            if (s[i] == 95)
                s[i] = '*';
            s_ori[s_ori_len++] = s[i];
        }
        else
        {
            s[i] = '-';
            res->gaps++;
            if (!ga)
            {
                ga = true;
                res->gap_open++;
            }
        }
    }

    if (has_must_include)
    {
        string s_ori_str(s_ori, len);
        if (!check_include(s_ori_str))
        {
            res->report = false;
            return;
        }
    }
    ga = false;
    for (int i = 0; i < len; i++)
    {
        // cout<<q_res[t][i]<<" ";
        if (res->q_res[i] != (size_t)(-1))
        {
            ga = false;
            q_seq[i] = q[res->q_res[i]] + 65;
        }
        else
        {
            q_seq[i] = '-';
            res->gaps++;
            if (!ga)
            {
                ga = true;
                res->gap_open++;
            }
        }
    }
    res->mismatch = 0;
    res->positive = 0;
    for (int i = 0; i < len; i++)
    {
        match[i]=' ';
        if (BLOSUM62[(q_seq[i] - 65) * 26 + (s[i] - 65)] > 0 && q_seq[i]!='-' && s[i]!='-')
        {
            res->positive++;
            match[i]='+';
        }
        if (q_seq[i] != s[i])
        {
            res->mismatch++;
        }
        else
        {
            match[i]=q_seq[i];
        }
            
    }
    res->n_identity = res->align_length - res->mismatch;
    res->p_identity = (1 - (double)res->mismatch / res->align_length) * 100;
    if (detailed_alignment)
    {
        res->q = q_seq;
        res->s = s;
        res->s_ori = s_ori;
        res->match = match;
    }
    
}

// void banded_smith_waterman(const char *q, const char *c, vector<uint32_t> &q_idxs, vector<uint32_t> &q_lens, vector<size_t> &diags, size_t c_len, size_t num_task, vector<SWResult> &res, ThreadPool *pool, vector<future<int>> &rs)
// {
// mu1.lock();

// vector<future<int>> rs;
// std::thread threads[num_task];

// mu1.lock();

// for (int i = 0; i < num_task; ++i)
// {
//     rs.emplace_back((*pool).enqueue([=, &res]
//                                  {
//         banded_smith_waterman_kernel(i,q,c,q_idxs[i],q_lens[i],diags[i],c_len,res);
//         return i; }));
//     // threads[i] = std::thread(banded_smith_waterman_kernel, i, q, c, q_idxs, q_lens, diags, c_len, ref(res));
//     // boost::asio::post(*pool, [=,&res](){
//     //     banded_smith_waterman_kernel(i, q, c, q_idxs, q_lens, diags, c_len, res);
//     // });
// }

// mu1.unlock();

// for (auto &r : rs)
// {
//     r.get();
// }

// for (auto &thread : threads)
// {
//     thread.join();
// }

// mu1.unlock();
// }