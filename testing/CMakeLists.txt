

CREATE_TEST_SOURCELIST(Tests MathTests.cpp
  TestGmres.cpp
  TestInexactNewton.cpp
)

ADD_EXECUTABLE(MathTests ${Tests})

ADD_TEST(TestGmres MathTests TestGmres)
ADD_TEST(TestInexactNewton MathTests TestInexactNewton)
