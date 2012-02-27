project(exafmm)

#FindMPI
option(USE_MPI "Use MPI" ON)
IF(USE_MPI)
    find_package(MPI REQUIRED)
    MESSAGE(STATUS "Enabling MPI")
    include_directories(${MPI_CXX_INCLUDE_PATH})
ENDIF(USE_MPI)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)
ADD_DEFINITIONS(-D_FILE_OFFSET_BITS=64)

option(USE_CPU "Use CPU kernels" OFF)
option(USE_GPU "Use CPU kernels" OFF)
option(USE_CPU_GPU "Use CPU kernels" OFF)
option(USE_SPHERICAL "Use Spherical Harmonics Expansions" OFF)
option(USE_CARTESIAN "Use Cartesian Expansions" OFF)

if(USE_CPU)
    set(USE_GPU OFF)
    set(USE_CPU_GPU OFF)
elseif(USE_GPU)
    set(USE_CPU OFF)
    set(USE_CPU_GPU OFF)
elseif(USE_CPU_GPU)
    set(USE_GPU OFF)
    set(USE_CPU OFF)
endif(USE_CPU)

if(USE_SPHERICAL)
    set(USE_CARTESIAN OFF)    
endif(USE_SPHERICAL)

if(USE_CARTESIAN)
    set(USE_SPHERICAL OFF)    
endif(USE_CARTESIAN)

if(USE_CPU)
    add_definitions(-Dcpu)
endif(USE_CPU)

if(USE_GPU)
    add_definitions(-Dgpu)
endif(USE_GPU)

if(USE_CPU_GPU)
    add_definitions(-Dcpugpu)
endif(USE_CPU_GPU)

if(USE_SPHERICAL)
    add_definitions(-DSpherical)
endif(USE_SPHERICAL)

if(USE_CARTESIAN)
    add_definitions(-DCartesian)
endif(USE_CARTESIAN)

option(USE_QUARK "Use Quark Tasking Library" OFF)
if(USE_QUARK)
    MESSAGE(STATUS "Enabling Quark")
    add_definitions(-DQUARK)
endif(USE_QUARK)

INCLUDE_DIRECTORIES(${CUDA_SDK_ROOT_DIR}/C/common/inc)

add_subdirectory(kernel)
add_subdirectory(fast)
add_subdirectory(external)
IF(BUILD_TESTING)
    add_subdirectory(unit_test)
ENDIF(BUILD_TESTING)
