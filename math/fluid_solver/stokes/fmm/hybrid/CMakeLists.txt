
#FindCUDA
IF(USE_CUDA)
  SET(CUDA_TARGETS_SOURCES kernel.cu)
  cuda_add_library(kernel STATIC ${CUDA_TARGETS_SOURCES})
  TARGET_LINK_LIBRARIES(kernel ${CUDA_LIBRARIES})
  add_library(hybridfmm hybridfmm.c expansions.c tree.c)
  target_link_libraries(hybridfmm kernel)
ENDIF(USE_CUDA)
# add_executable(main main.cpp)
# target_link_libraries(main hybridfmm)

# FIND_PACKAGE(Qt4 REQUIRED)
# SET(QT_USE_QTOPENGL TRUE)
#
IF(USE_QT_GUI)
    QT4_WRAP_CPP(renderer_HEADERS_MOC tree_renderer.hpp)
    add_library(octreerenderer ${renderer_HEADERS_MOC})
    target_link_libraries(octreerenderer ${QT_LIBRARIES} GL GLU)
ENDIF(USE_QT_GUI)


