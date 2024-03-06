cmake_minimum_required(VERSION 3.5)

set(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda)
project(BLASTP LANGUAGES CXX C CUDA)

# SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_BUILD_TYPE DEBUG)

option(GLF_GPU_SW "Use GLF-GPU Smith-waterman" ON)
if(GLF_GPU_SW)
    add_definitions(-DGLF_GPU_SW)
endif()

option(USE_GPU_SW "Use GPU Smith-waterman" OFF)
if(USE_GPU_SW)
    add_definitions(-DUSE_GPU_SW)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++17 -Xcompiler -pthread")

if (NOT ("${CMAKE_SIZEOF_VOID_P}" STREQUAL "8"))
    message(SEND_ERROR "require 64 bit system")
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "DEBUG MODE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O0 -lineinfo -g")
else()
    message(STATUS "RELEASE MODE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")
endif()

find_package(OpenMP)

message(STATUS "Use CPU local alignment")
set(SOURCES_UTIL
	    src/banded_smith_cpu.cpp
        src/util.cpp
        src/qit.cpp
        src/output.cpp
)

set(SOURCES_SEARCH
    src/hash_table.cu
    src/blastp.cu
    src/query_group.cu
)

add_library(util STATIC ${SOURCES_UTIL})
target_include_directories(util  PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(util PRIVATE ${OpenMP_CXX_LIBRARIES})
target_include_directories (util
	PRIVATE
		 ${OpenMP_CXX_INCLUDE_DIRS}
)
target_compile_options (util
	PRIVATE
		${OpenMP_CXX_FLAGS}
)

# gpu_sw
add_subdirectory(./gpu-sw)

add_executable(query src/main.cpp ${SOURCES_SEARCH})

# find_library(GPU_SW_LIBRARY
#     PATHS ${CMAKE_BINARY_DIR}/gpu-sw
# )
# target_link_libraries(query PRIVATE util ${GPU_SW_LIBRARY})
target_link_libraries(query PRIVATE util sw-lib)
# target_link_libraries(query util)

add_executable(createDB src/createDB.cpp)

target_include_directories(query
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/gpu-sw
)

target_include_directories(createDB
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)

