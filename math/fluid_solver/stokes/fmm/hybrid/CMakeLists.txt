
#FindCUDA
FIND_PACKAGE(CUDA REQUIRED)
SET(CUDA_NVCC_FLAGS --compiler-bindir /opt/gcc-4.4.6/bin/ -gencode arch=compute_10,code=sm_10)


SET(CUDA_TARGETS_SOURCES kernel.cu)
cuda_add_library(kernel STATIC ${CUDA_TARGETS_SOURCES})
TARGET_LINK_LIBRARIES(kernel ${CUDA_LIBRARIES})

add_library(hybridfmm hybridfmm.c expansions.c tree.c)
target_link_libraries(hybridfmm kernel)

add_executable(main main.cpp)
target_link_libraries(main hybridfmm)
