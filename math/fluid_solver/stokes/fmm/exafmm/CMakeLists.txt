project(exafmm)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)
ADD_DEFINITIONS(-DSpherical -Dcpu -D_FILE_OFFSET_BITS=64)

INCLUDE_DIRECTORIES(${CUDA_SDK_ROOT_DIR}/C/common/inc)

# add_custom_target(
#   doc-prerequisites
#   ALL
#   COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs/
#   COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs/html/
#   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
#                                    ${CMAKE_CURRENT_BINARY_DIR}/docs
#   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
# )
# 
# add_custom_target(doc ALL
#   COMMAND doxygen
#   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs)
# 
# add_dependencies(doc doc-prerequisites)

add_subdirectory(kernel)
add_subdirectory(fast)
add_subdirectory(unit_test)
add_subdirectory(lib)