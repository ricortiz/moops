
#FindCUDA
FIND_PACKAGE(CUDA REQUIRED)
SET(CUDA_NVCC_FLAGS --compiler-bindir /usr/bin/gcc-4.4 -gencode arch=compute_20,code=sm_20)


SET(CUDA_TARGETS_SOURCES kernel.cu)
cuda_add_library(kernel STATIC ${CUDA_TARGETS_SOURCES})
TARGET_LINK_LIBRARIES(kernel ${CUDA_LIBRARIES})

add_library(hybridfmm hybridfmm.c expansions.c tree.c)
target_link_libraries(hybridfmm kernel)

add_executable(main main.cpp)
target_link_libraries(main hybridfmm)

FIND_PACKAGE(Qt4 REQUIRED)
SET(QT_USE_QTOPENGL TRUE)

INCLUDE(${QT_USE_FILE})
QT4_WRAP_CPP(renderer_HEADERS_MOC tree_renderer.hpp)

add_library(otgui ${renderer_HEADERS_MOC})
target_link_libraries(otgui ${QT_LIBRARIES} GL GLU hybridfmm)
