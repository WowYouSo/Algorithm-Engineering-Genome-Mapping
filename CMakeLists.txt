cmake_minimum_required(VERSION 4.1)
project(Genome_Mapping)

set(CMAKE_CXX_STANDARD 20)

if(APPLE)
    set(OpenMP_C_FLAGS "-Xclang -fopenmp")
    set(OpenMP_CXX_FLAGS "-Xclang -fopenmp")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    include_directories("/opt/homebrew/opt/libomp/include")
endif()

find_package(OpenMP REQUIRED)

add_executable(Genome_Mapping main.cpp)
target_link_libraries(Genome_Mapping OpenMP::OpenMP_CXX)
