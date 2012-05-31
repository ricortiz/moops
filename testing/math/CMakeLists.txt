
SET(fluid_solvers cpu_stokes_solver.cpp images.cpp)
SET(linear_solvers generalized_minimal_residual_method.cpp)
SET(nonlinear_solvers inexact_newton.cpp)
SET(ode_solvers backward_euler.cpp forward_euler.cpp explicit_sdc.cpp semi_implicit_sdc.cpp)

if(USE_EXAFMM)
    SET(fluid_solvers ${fluid_solvers} exafmm_stokes_solver.cpp)
    ADD_DEFINITIONS(-DSTOKES)
    ADD_DEFINITIONS(-DSPHERICAL)
endif()

make_tests("${fluid_solvers}" "TestFluidSolvers")
make_tests("${linear_solvers}" "TestLinearSolvers")
make_tests("${nonlinear_solvers}" "TestNonLinearSolvers")
make_tests("${ode_solvers}" "TestOdeSolvers")

TARGET_LINK_LIBRARIES(TestFluidSolvers ${APP_LIBS})
if(USE_EXAFMM)
    TARGET_LINK_LIBRARIES(TestFluidSolvers Kernels)
    IF(USE_QUARK)
        TARGET_LINK_LIBRARIES(TestFluidSolvers quark)
    ENDIF(USE_QUARK)
endif()

add_executable(conovergence_rate conovergence_rate.cpp)
target_link_libraries(conovergence_rate mgl ${APP_LIBS})


