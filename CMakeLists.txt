cmake_minimum_required(VERSION 2.8)

project(MOOPS)

SET(APP_LIBS "" CACHE INTERNAL "")
#FindVTK
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
SET(APP_LIBS ${APP_LIBS} QVTK vtkRendering vtkGraphics vtkIO vtkCommon)
# FindQt4
FIND_PACKAGE(Qt4 REQUIRED)
SET(QT_USE_QTOPENGL TRUE)
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

SET(USE_CUDA_FLUID_SOLVER OFF CACHE BOOL "Use the cuda-based Stokes solver.")

include_directories(${PROJECT_SOURCE_DIR})

add_subdirectory(app)
IF(USE_CUDA_FLUID_SOLVER)
  add_subdirectory(math/fluid_solver/stokes/nbody_gpu)
ENDIF(USE_CUDA_FLUID_SOLVER)

add_subdirectory(gui)