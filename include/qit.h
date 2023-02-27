#pragma once

#include "cmath"
#include "params.h"
#include <string>
#include <iostream>

using namespace std;

class QIT
{

public:
    uint8_t *index_size; // LENGTH
    int *q_num;          // LENGTH * WIDTH
    int *q_idx;          // LENGTH * WIDTH

    QIT();
    ~QIT();

    inline size_t BitsToBytes(size_t bits) { return (bits - 1) / 8 + 1; }

    inline int get_mers_num(string mers)
    {
        int sum = 0;
        for (int i = 0; i < seed_length; i++)
        {
            sum += mers[i] << (5 * i);
        }
        return sum;
    }

    int get_index_size(string mers);

    void add(string mers, int q_num, int q_idx);

    int find_empty_dummy(int loc);
};
