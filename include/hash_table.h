#pragma once

#include "params.h"
#include "util.h"

struct KeyValue
{
    uint32_t key;
    uint32_t value;
};

struct Task
{
    uint32_t key;
    uint16_t value;
    uint16_t q_id;
};

const uint32_t kEmpty = 0xffffffff;

KeyValue *create_hashtable(size_t size);

void destroy_hashtable(KeyValue *hashtable);
