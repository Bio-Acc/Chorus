#include "qit.h"
#include "hash_table.h"
#include "smith.h"
#include "util.h"
#include "query_group.h"
#include "output.h"

struct TimeProfile 
{
    double mem_time = 0;
    double gpu_time = 0;
    double cpu_time = 0;
    double name_time = 0;
};

void blastp(string argv_query, vector<string> argv_dbs, string argv_out);