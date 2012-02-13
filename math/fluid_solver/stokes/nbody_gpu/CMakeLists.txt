
#FindCUDA
IF(USE_CUDA_FLUID_SOLVER)
  FIND_PACKAGE(CUDA REQUIRED)
  SET(CUDA_TARGETS_SOURCES kernels/ComputeImages.cu kernels/ComputeImagesKernel.cu kernels/ComputeStokeslets.cu kernels/ComputeStokesletsKernel.cu kernels/ComputeStokesletsNearField.cu kernels/ComputeStokesletsNearFieldKernel.cu kernels/ComputeImagesNearField.cu kernels/ComputeImagesNearFieldKernel.cu)
  SET(CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
  cuda_add_library(cudaStokesSolver STATIC ${CUDA_TARGETS_SOURCES})
  TARGET_LINK_LIBRARIES(cudaStokesSolver ${CUDA_LIBRARIES})
ENDIF(USE_CUDA_FLUID_SOLVER)
