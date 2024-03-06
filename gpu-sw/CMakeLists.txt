add_library(sw-lib STATIC)
set_property(TARGET sw-lib PROPERTY CUDA_ARCHITECTURES "${GPU_ARCHS}")
target_sources(sw-lib PRIVATE
	"sw.cu"
	"cuda-utils.cc"
	"mats.cc"
	"seqs.cc"
	"sw-lib.cc"
	)
target_link_libraries(sw-lib PUBLIC "${CUDA_LIBRARIES}")
target_include_directories(sw-lib PUBLIC "${CUDA_INCLUDE_DIRS}" ".")