
add_library(cpuKernel cpuEwaldLaplace.cxx cpuP2P.cxx cpuSphericalLaplace.cxx cpuVanDerWaals.cxx)

SET(CUDA_TARGETS_SOURCES gpuSphericalLaplace.cu gpuVanDerWaals.cu cpuP2P.cxx)
cuda_add_library(gpuKernel STATIC ${CUDA_TARGETS_SOURCES})
TARGET_LINK_LIBRARIES(gpuKernel ${CUDA_LIBRARIES})
