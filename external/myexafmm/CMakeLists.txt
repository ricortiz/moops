project(exafmm)
cmake_minimum_required(VERSION 2.8)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)
ADD_DEFINITIONS(-D_FILE_OFFSET_BITS=64)

FIND_PACKAGE(VTK REQUIRED)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
    SET(APP_LIBS ${APP_LIBS} vtkIO QVTK vtkRendering vtkGraphics vtkCommon)
    MESSAGE(STATUS "Enabling VTK")
ENDIF(VTK_FOUND)

#FindMPI
option(USE_MPI "Use MPI" ON)
IF(USE_MPI)
    FIND_PACKAGE(MPI REQUIRED)
    set(CMAKE_CXX_COMPILER /usr/local/bin/mpic++)
    MESSAGE(STATUS "Enabling MPI")
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
ENDIF()

option(USE_CPU "Use CPU kernels" OFF)
option(USE_GPU "Use CPU kernels" OFF)
option(USE_SPHERICAL "Use Spherical Harmonics Expansions" ON)
option(USE_CARTESIAN "Use Cartesian Expansions" OFF)

if(USE_CPU)
    add_definitions(-Dcpu)
endif(USE_CPU)

if(USE_GPU)
  ADD_DEFINITIONS(-Dgpu)
  FIND_PACKAGE(CUDA REQUIRED)
  SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode arch=compute_20,code=sm_20 -use_fast_math)
  IF(NOT CMAKE_COMPILER_IS_GNUCXX)
    SET(CUDA_PROPAGATE_HOST_FLAGS OFF)
  ENDIF()
  MESSAGE(STATUS "Enabling CUDA ${CUDA_VERSION}")
  IF(${CUDA_VERSION} MATCHES "4.1")
    ADD_DEFINITIONS(-DCUDA_4_1)
  ENDIF()
endif(USE_GPU)

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

INCLUDE(CTest)
IF(BUILD_TESTING)
    add_subdirectory(unit_test)
ENDIF(BUILD_TESTING)

add_subdirectory(kernel)
add_subdirectory(fast)
add_subdirectory(external)

