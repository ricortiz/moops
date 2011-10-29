
set(CUDA_TARGETS_SOURCES kernels/ComputeImages.cu kernels/ComputeImagesKernel.cu kernels/ComputeStokeslets.cu kernels/ComputeStokesletsKernel.cu kernels/ComputeStokesletsNearField.cu kernels/ComputeStokesletsNearFieldKernel.cu kernels/ComputeImagesNearField.cu kernels/ComputeImagesNearFieldKernel.cu)
#FindCUDA
find_package(CUDA)
IF(CUDA_FOUND)
    set(CUDA_NVCC_FLAGS --compiler-bindir /opt/gcc-4.4.6/bin -gencode arch=compute_20,code=sm_20)
    cuda_add_library(cudaStokesSolver STATIC ${CUDA_TARGETS_SOURCES})
    target_link_libraries(cudaStokesSolver ${CUDA_LIBRARIES})
ENDIF(CUDA_FOUND)
