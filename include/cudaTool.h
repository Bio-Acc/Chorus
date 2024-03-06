#include <iostream>
#include <cuda_runtime.h>
/**
 * @brief get the avaialbe VRAM of GPUs
 * @param defaultIndexOfGPU the index of GPUs
**/
double getVRAM(int defaultIndexOfGPU = 0);