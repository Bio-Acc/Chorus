#include "cudaTool.h"

#define CUDA_CALL(F)                                                          \
    if ((F) != cudaSuccess)                                                   \
    {                                                                         \
        printf("Error %s at %s:%d\n", cudaGetErrorString(cudaGetLastError()), \
               __FILE__, __LINE__);                                           \
        exit(-1);                                                             \
    }
double getVRAM(int defaultIndexOfGPU)
{
    int deviceCount;
    CUDA_CALL(cudaGetDeviceCount(&deviceCount));
    if (defaultIndexOfGPU >= deviceCount) {
        std::cerr << "No CUDA-enabled devices found." << std::endl;
        return 1;
    }
    CUDA_CALL(cudaSetDevice(defaultIndexOfGPU));
    size_t freeMem, totalMem;
    CUDA_CALL(cudaMemGetInfo(&freeMem, &totalMem));
    double result = (double)freeMem / (1024 * 1024 * 1024) / 3;
    return max(1.0, result); 
}