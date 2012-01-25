
#FindCUDA
FIND_PACKAGE(CUDA REQUIRED)
SET(CUDA_TARGETS_SOURCES kernel.cu)
SET(CUDA_NVCC_FLAGS --compiler-bindir /opt/bin -gencode arch=compute_20,code=sm_20)
cuda_add_library(kernel STATIC kernel.cu)
TARGET_LINK_LIBRARIES(kernel ${CUDA_LIBRARIES})



add_library(hybridfmm driver.c expansions.c tree.c)
target_link_libraries(hybridfmm kernel)
