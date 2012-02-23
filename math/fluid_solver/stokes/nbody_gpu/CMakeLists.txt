
#FindCUDA
IF(USE_CUDA)
  SET(CUDA_TARGETS_SOURCES kernels/ComputeStokesletsKernel.cu kernels/ComputeStokeslets.cu kernels/ComputeStokesletsNearFieldKernel.cu)
  cuda_add_library(cudaStokesSolver STATIC ${CUDA_TARGETS_SOURCES})
  TARGET_LINK_LIBRARIES(cudaStokesSolver ${CUDA_LIBRARIES})
ENDIF(USE_CUDA)
