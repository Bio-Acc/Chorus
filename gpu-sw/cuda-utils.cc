#include "cuda-utils.hh"

#include <cstdio>
#include <iostream>

void eccl::dump_device_info(int device) {
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, device);
	printf("%s\n", prop.name);
	printf("Major revision number:         %d\n", prop.major);
	printf("Minor revision number:         %d\n", prop.minor);
	printf("Total global memory:           %zu", prop.totalGlobalMem);
	printf(" bytes\n");
	printf("Number of multiprocessors:     %d\n", prop.multiProcessorCount);
	printf("Total amount of shared memory per block: %zu\n",prop.sharedMemPerBlock);
	printf("Total registers per block:     %d\n", prop.regsPerBlock);
	printf("Warp size:                     %d\n", prop.warpSize);
	printf("Maximum memory pitch:          %zu\n", prop.memPitch);
	printf("Total amount of constant memory:         %zu\n",   prop.totalConstMem);
}

namespace eccl {
void operator,(cudaError_t error, eccl::check_cuda checker) {
	if(!error)
		return;
	std::cerr<<"error: "<<checker._msg<<": "<<cudaGetErrorName(error)<<": "<<cudaGetErrorString(error)<<"\n";
	std::exit(EXIT_FAILURE);
}
}
