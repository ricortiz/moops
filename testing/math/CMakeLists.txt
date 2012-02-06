
set(fluid_solvers cpu_stokes_solver fmm_stokes_solver hybrid_fmm_stokes_solver kernel)
set(linear_solvers generalized_minimal_residual_method)
set(nonlinear_solvers inexact_newton)
set(ode_solvers backward_euler forward_euler explicit_sdc semi_implicit_sdc)
set(target_links hybridfmm)
FIND_PACKAGE(Qt4 REQUIRED)
SET(QT_USE_QTOPENGL TRUE)

INCLUDE(${QT_USE_FILE})
make_tests("${fluid_solvers}" "TestFluidSolvers" "otgui")
make_tests("${linear_solvers}" "TestLinearSolvers" "")
make_tests("${nonlinear_solvers}" "TestNonLinearSolvers" "")
make_tests("${ode_solvers}" "TestOdeSolvers" "")
