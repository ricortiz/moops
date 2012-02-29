

add_library(cpuKernel cpuP2P.cxx cpuSphericalStokes.cxx)

if(USE_CUDA)
    SET(CUDA_TARGETS_SOURCES gpuSphericalLaplace.cu gpuSphericalStokes.cu cpuP2P.cxx)
    cuda_add_library(gpuKernel STATIC ${CUDA_TARGETS_SOURCES})
    TARGET_LINK_LIBRARIES(gpuKernel ${CUDA_LIBRARIES})
endif(USE_CUDA)
