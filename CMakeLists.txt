cmake_minimum_required(VERSION 2.8)

project(MOOPS)

SET(APP_LIBS "" CACHE INTERNAL "")
SET(USE_CUDA_FLUID_SOLVER OFF CACHE BOOL "Use the cuda-based Stokes solver.")
SET(USE_QT_GUI OFF CACHE BOOL "Build the QT4 GUI framework.")

#FindVTK
FIND_PACKAGE(VTK)
IF(VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
  SET(APP_LIBS ${APP_LIBS} vtkIO)
ENDIF(VTK_FOUND)

# FindQt4
if(USE_QT_GUI)
    FIND_PACKAGE(Qt4 REQUIRED)
    SET(QT_USE_QTOPENGL TRUE)
    INCLUDE(${QT_USE_FILE})
    ADD_DEFINITIONS(${QT_DEFINITIONS} -DUSE_QT_GUI)
    SET(APP_LIBS ${APP_LIBS} gui QVTK vtkRendering vtkGraphics vtkCommon)
    ADD_SUBDIRECTORY(gui)
endif(USE_QT_GUI)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/testing)

include(CheckCXXCompilerFlag)
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++ -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -fexceptions -fno-check-new -fno-common -fstrict-aliasing")
    set(CMAKE_CXX_FLAGS_DEBUG "-g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-g0 -O2")
    check_cxx_compiler_flag("-Wno-variadic-macros" COMPILER_SUPPORT_WNOVARIADICMACRO)
    if(COMPILER_SUPPORT_WNOVARIADICMACRO)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-variadic-macros")
    endif()

    check_cxx_compiler_flag("-Wextra" COMPILER_SUPPORT_WEXTRA)
    if(COMPILER_SUPPORT_WEXTRA)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    endif()

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")

    check_cxx_compiler_flag("-fopenmp" COMPILER_SUPPORT_OPENMP)
    if(COMPILER_SUPPORT_OPENMP)
        option(ENABLE_OPENMP_GLOBALY "Enable/Disable OpenMP" OFF)
        if(ENABLE_OPENMP_GLOBALY)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            message(STATUS "Enabling OpenMP")
        endif()
    endif()
endif(CMAKE_COMPILER_IS_GNUCXX)

IF(USE_CUDA_FLUID_SOLVER)
  ADD_SUBDIRECTORY(math/fluid_solver/stokes/nbody_gpu/)
  ADD_DEFINITIONS(-DUSE_CUDA_FLUID_SOLVER)
  SET(APP_LIBS ${APP_LIBS} cudaStokesSolver)
ENDIF(USE_CUDA_FLUID_SOLVER)

INCLUDE(CTest)
IF(BUILD_TESTING)
    ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)


option(ENABLE_EXAMPLES "Enable examples" ON)
if(ENABLE_EXAMPLES)
        add_subdirectory(examples)
endif(ENABLE_EXAMPLES)

ADD_SUBDIRECTORY(math/fluid_solver/stokes/fmm/hybrid)
# 