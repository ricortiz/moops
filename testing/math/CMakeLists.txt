
set(fluid_solvers cpu_stokes_solver fmm_stokes_solver hybrid_fmm_stokes_solver)
set(linear_solvers generalized_minimal_residual_method)
set(nonlinear_solvers inexact_newton)
set(ode_solvers backward_euler forward_euler explicit_sdc semi_implicit_sdc)


make_tests("${fluid_solvers}" "TestFluidSolvers")
make_tests("${linear_solvers}" "TestLinearSolvers")
make_tests("${nonlinear_solvers}" "TestNonLinearSolvers")
make_tests("${ode_solvers}" "TestOdeSolvers")

target_link_libraries(TestFluidSolvers ${APP_LIBS})
IF(USE_QT_GUI)
    target_link_libraries(TestFluidSolvers octreerenderer)
ENDIF(USE_QT_GUI)
