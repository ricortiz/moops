project(QUARK)

find_package(Threads REQUIRED)

add_library(quark STATIC icl_list.c icl_hash.c quark.c quarkos.c)
target_link_libraries(quark pthread)

file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/include)

install(FILES bsd_queue.h bsd_tree.h icl_hash.h icl_list.h quark.h quarkwinthread.h DESTINATION ${PROJECT_SOURCE_DIR}/include)

include_directories(${PROJECT_SOURCE_DIR}/include)

