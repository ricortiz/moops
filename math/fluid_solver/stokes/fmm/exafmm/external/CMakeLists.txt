set(ext_libs quark)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ext_include)

FOREACH(dir ${ext_libs})
    file(GLOB_RECURSE files ./${dir} *.h *.hpp)
    file(COPY ${files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ext_include)
ENDFOREACH(dir)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/ext_include)
if(USE_QUARK)
    add_library(quark quark/icl_list.c quark/icl_hash.c quark/quark.c quark/quarkos.c)
    target_link_libraries(quark pthread)
endif(USE_QUARK)
