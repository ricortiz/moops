cmake_minimum_required(VERSION 2.8)

project(oopsi)

FIND_PACKAGE(VTK REQUIRED)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
    SET(APP_LIBS ${APP_LIBS} vtkIO)
ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR "VTK library not found.")
ENDIF(VTK_FOUND)



FIND_PACKAGE(Qt4 REQUIRED)
SET(QT_USE_QTOPENGL TRUE)
# 
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

SET(USE_CUDA_FLUID_SOLVER OFF CACHE BOOL "Use the cuda-based Stokes solver.")

include_directories(${PROJECT_SOURCE_DIR})

QT4_WRAP_CPP(renderer_HEADERS_MOC io/particle_renderer.hpp)
add_executable(renderer renderer.cpp ${renderer_HEADERS_MOC})
target_link_libraries(renderer ${QT_LIBRARIES} )

add_subdirectory(app)
add_subdirectory(math)
add_subdirectory(io)
