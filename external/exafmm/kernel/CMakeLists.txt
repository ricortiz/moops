

add_library(cpuKernel cpuP2P.cxx cpuSphericalLaplace.cxx)

if(USE_GPU)
    SET(CUDA_TARGETS_SOURCES gpuSphericalLaplace.cu cpuP2P.cxx)
    cuda_add_library(gpuKernel STATIC ${CUDA_TARGETS_SOURCES})
    TARGET_LINK_LIBRARIES(gpuKernel ${CUDA_LIBRARIES})
endif(USE_GPU)
