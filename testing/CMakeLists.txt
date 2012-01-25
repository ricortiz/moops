

# CREATE_TEST_SOURCELIST(Tests MathTests.cpp
#   TestGmres.cpp
#   TestInexactNewton.cpp
#   EXTRA_INCLUDE TestUtilities.h
# )
# 
# ADD_EXECUTABLE(MathTests ${Tests})
# 
# SET (TestsToRun ${Tests})
# REMOVE (TestsToRun MathTests.cpp)
# 
# FOREACH (test ${TestsToRun})
#     GET_FILENAME_COMPONENT(TName ${test} NAME_WE)
#     ADD_TEST(${TName} ${CXX_TEST_PATH}/MathTests ${TName})
# ENDFOREACH (test)
