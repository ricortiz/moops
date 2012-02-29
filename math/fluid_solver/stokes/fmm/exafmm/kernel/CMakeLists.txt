

add_library(cpuKernel cpuEwaldLaplace.cxx cpuP2P.cxx cpuSphericalLaplace.cxx cpuSphericalStokes.cxx cpuVanDerWaals.cxx)

if(USE_CUDA)
SET(CUDA_TARGETS_SOURCES gpuSphericalLaplace.cu gpuSphericalStokes.cu gpuVanDerWaals.cu cpuP2P.cxx)
cuda_add_library(gpuKernel STATIC ${CUDA_TARGETS_SOURCES})
TARGET_LINK_LIBRARIES(gpuKernel ${CUDA_LIBRARIES})
endif(USE_CUDA)
