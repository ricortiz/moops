
IF(USE_QUARK)
    INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/math/fluid_solver/stokes/fmm/exafmm/external/include)
    add_definitions(-DQUARK)
    set(APP_LIBS ${APP_LIBS} quark)
ENDIF(USE_QUARK)
add_definitions(-DSpherical)
if(USE_CUDA)
    add_definitions(-Dcpugpu)
else(USE_CUDA)
    add_definitions(-Dcpu)
endif(USE_CUDA)

set(fluid_solvers cpu_stokes_solver.cpp images.cpp exafmm_stokes_solver.cpp)
set(linear_solvers generalized_minimal_residual_method.cpp)
set(nonlinear_solvers inexact_newton.cpp)
set(ode_solvers backward_euler.cpp forward_euler.cpp explicit_sdc.cpp semi_implicit_sdc.cpp)

set(spherical spherical)

make_tests("${spherical}" "TestSphericalHarmonics")
make_tests("${fluid_solvers}" "TestFluidSolvers")
make_tests("${linear_solvers}" "TestLinearSolvers")
make_tests("${nonlinear_solvers}" "TestNonLinearSolvers")
make_tests("${ode_solvers}" "TestOdeSolvers")

target_link_libraries(TestFluidSolvers ${APP_LIBS} cpuKernel)

IF(USE_QT_GUI)
    target_link_libraries(TestFluidSolvers octreerenderer)
ENDIF(USE_QT_GUI)
