
add_library(quark quark/icl_list.c quark/icl_hash.c quark/quark.c quark/quarkos.c)
target_link_libraries(quark pthread)
