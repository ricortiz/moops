
IF(USE_QUARK)
    INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/math/fluid_solver/stokes/fmm/exafmm/external/include)
    ADD_DEFINITIONS(-DQUARK)
    SET(APP_LIBS ${APP_LIBS} quark)
ENDIF(USE_QUARK)
ADD_DEFINITIONS(-DSpherical)


SET(fluid_solvers cpu_stokes_solver.cpp images.cpp stokesParallel.cxx)
SET(linear_solvers generalized_minimal_residual_method.cpp)
SET(nonlinear_solvers inexact_newton.cpp)
SET(ode_solvers backward_euler.cpp forward_euler.cpp explicit_sdc.cpp semi_implicit_sdc.cpp)

SET(spherical spherical)

make_tests("${spherical}" "TestSphericalHarmonics")
make_tests("${fluid_solvers}" "TestFluidSolvers")
make_tests("${linear_solvers}" "TestLinearSolvers")
make_tests("${nonlinear_solvers}" "TestNonLinearSolvers")
make_tests("${ode_solvers}" "TestOdeSolvers")

TARGET_LINK_LIBRARIES(TestFluidSolvers ${APP_LIBS})

IF(USE_QT_GUI)
    TARGET_LINK_LIBRARIES(TestFluidSolvers octreerenderer)
ENDIF(USE_QT_GUI)
