cmake_minimum_required(VERSION 3.5)

if (POLICY CMP0104)
	cmake_policy(SET CMP0104 OLD)
endif()

set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
project(BLASTP LANGUAGES CXX CUDA)

SET(CMAKE_BUILD_TYPE "Release")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread")
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

include_directories(
    ${PROJECT_SOURCE_DIR}/include
)

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
target_link_libraries(util PRIVATE ${OpenMP_CXX_LIBRARIES})
target_include_directories (util
	PRIVATE ${OpenMP_CXX_INCLUDE_DIRS}
)
target_compile_options (util
	PRIVATE ${OpenMP_CXX_FLAGS}
)

add_executable(query src/main.cpp ${SOURCES_SEARCH})
target_link_libraries(query util)

add_executable(createDB src/createDB.cpp)

