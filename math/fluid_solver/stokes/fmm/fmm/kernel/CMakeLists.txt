
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/math/fluid_solver/stokes/fmm/fmm/include)
INCLUDE_DIRECTORIES(${CUDA_SDK_ROOT_DIR}/C/common/inc)
add_definitions(-DSpherical)
if(USE_CUDA)
    SET(CUDA_TARGETS_SOURCES gpuSphericalStokes.cu cpuSphericalStokes.cpp cpuP2P.cpp)
    cuda_add_library(hybridFMMStokesSolver STATIC ${CUDA_TARGETS_SOURCES})
    TARGET_LINK_LIBRARIES(hybridFMMStokesSolver ${CUDA_LIBRARIES})
endif(USE_CUDA)