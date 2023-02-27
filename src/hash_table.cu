#include "stdio.h"
#include "stdint.h"
#include "vector"
#include "hash_table.h"

// Create a hash table. For linear probing, this is just an array of KeyValues
KeyValue *create_hashtable(size_t size)
{
    cudaError_t cudaStatus;
    // Allocate memory
    KeyValue *hashtable;
    cudaStatus = cudaMalloc(&hashtable, sizeof(KeyValue) * size);
    if (cudaStatus != cudaSuccess)
    {
        fprintf(stderr, "create_hashtable cudaMalloc err %s\n", cudaGetErrorString(cudaStatus));
        return nullptr;
    }

    // Initialize hash table to empty
    static_assert(kEmpty == 0xffffffff, "memset expected kEmpty=0xffffffff");
    cudaStatus = cudaMemset(hashtable, 0xff, sizeof(KeyValue) * size);
    if (cudaStatus != cudaSuccess)
    {
        fprintf(stderr, "create_hashtable cudaMemset err %d\n", cudaStatus);
    }

    return hashtable;
}

// Iterate over every item in the hashtable; return non-empty key/values
// __global__ void gpu_iterate_hashtable(KeyValue *pHashTable, KeyValue *kvs, uint32_t *kvs_size)
// {
//     unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;

//     if (threadid < kHashTableCapacity_dev)
//     {
//         if (pHashTable[threadid].key != kEmpty)
//         {
//             uint32_t value = pHashTable[threadid].value;
//             if (value != kEmpty and value >= MIN_DIAG_HIT)
//             {
//                 uint32_t size = atomicAdd(kvs_size, 1);
//                 kvs[size] = pHashTable[threadid];
//             }
//         }
//     }
// }

// std::vector<KeyValue> iterate_hashtable(KeyValue *pHashTable)
// {
//     cudaError_t cudaStatus;
//     uint32_t *device_num_kvs;
//     cudaStatus = cudaMalloc(&device_num_kvs, sizeof(uint32_t));
//     if (cudaStatus != cudaSuccess)
//     {
//         fprintf(stderr, "iterate_hashtable cudaMalloc err %d\n", cudaStatus);
//     }
//     cudaStatus = cudaMemset(device_num_kvs, 0, sizeof(uint32_t));
//     if (cudaStatus != cudaSuccess)
//     {
//         fprintf(stderr, "iterate_hashtable cudaMemset err %d\n", cudaStatus);
//     }

//     KeyValue *device_kvs;
//     cudaStatus = cudaMalloc(&device_kvs, sizeof(KeyValue) * (kHashTableCapacity_host / 2));
//     if (cudaStatus != cudaSuccess)
//     {
//         fprintf(stderr, "iterate_hashtable cudaMalloc err %d\n", cudaStatus);
//     }

//     cudaStatus = cudaMemcpyToSymbol(kHashTableCapacity_dev, &kHashTableCapacity_host, sizeof(uint32_t));
//     if (cudaStatus != cudaSuccess)
//     {
//         fprintf(stderr, "blastp cudaMemcpyToSymbol err %d\n", cudaStatus);
//     }

//     int mingridsize;
//     int threadblocksize;
//     cudaStatus = cudaOccupancyMaxPotentialBlockSize(&mingridsize, &threadblocksize, gpu_iterate_hashtable, 0, 0);
//     if (cudaStatus != cudaSuccess)
//     {
//         fprintf(stderr, "iterate_hashtable cudaOccupancyMaxPotentialBlockSize err %d\n", cudaStatus);
//     }

//     int gridsize = (kHashTableCapacity_host + threadblocksize - 1) / threadblocksize;
//     gpu_iterate_hashtable<<<gridsize, threadblocksize>>>(pHashTable, device_kvs, device_num_kvs);

//     uint32_t num_kvs;
//     cudaMemcpy(&num_kvs, device_num_kvs, sizeof(uint32_t), cudaMemcpyDeviceToHost);

//     std::vector<KeyValue> kvs;
//     kvs.resize(num_kvs);

//     cudaMemcpy(kvs.data(), device_kvs, sizeof(KeyValue) * num_kvs, cudaMemcpyDeviceToHost);

//     cudaFree(device_kvs);
//     cudaFree(device_num_kvs);

//     return kvs;
// }

// Free the memory of the hashtable
void destroy_hashtable(KeyValue *pHashTable)
{
    cudaFree(pHashTable);
}
