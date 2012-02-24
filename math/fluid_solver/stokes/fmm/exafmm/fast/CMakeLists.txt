add_definitions(-DSpherical -DMANY)
add_executable(main main.cxx)
target_link_libraries(main quark)
add_subdirectory(testing)