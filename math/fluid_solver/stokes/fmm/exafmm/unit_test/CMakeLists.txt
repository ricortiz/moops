
set(LIBS_OBJ mpi mpi_cxx vtkRendering vtkGraphics vtkFiltering vtkViews vtkCommon vtkWidgets vtkIO quark)

add_executable(construct construct.cxx)
target_link_libraries(construct ${LIBS_OBJ} cpuKernel)

add_executable(kernel kernel.cxx)
target_link_libraries(kernel cpuKernel)

add_executable(ewald_direct ewald_direct.cxx)
target_link_libraries(ewald_direct cpuKernel)

add_executable(serialrun serialrun.cxx)
target_link_libraries(serialrun ${LIBS_OBJ} cpuKernel)

add_executable(unsort unsort.cxx)
target_link_libraries(unsort ${LIBS_OBJ} cpuKernel)

add_executable(ijserialrun ijserialrun.cxx)
target_link_libraries(ijserialrun ${LIBS_OBJ} cpuKernel quark)

add_executable(Nserial Nserial.cxx)
target_link_libraries(Nserial ${LIBS_OBJ} cpuKernel)

# add_definitions(-DTOPDOWN)
add_executable(direct_gpu direct_gpu.cxx)
target_link_libraries(direct_gpu gpuKernel)

add_executable(mpi_test mpi.cxx)
target_link_libraries(mpi_test ${LIBS_OBJ})

add_executable(check_gpus check_gpus.cxx)
target_link_libraries(check_gpus ${LIBS_OBJ} cpuKernel)

add_executable(shift shift.cxx)
target_link_libraries(shift ${LIBS_OBJ} cpuKernel)

add_executable(nth_element nth_element.cxx)
target_link_libraries(nth_element ${LIBS_OBJ} cpuKernel)

add_executable(bisection bisection.cxx)
target_link_libraries(bisection ${LIBS_OBJ} cpuKernel)

add_executable(let let.cxx)
target_link_libraries(let ${LIBS_OBJ} cpuKernel)

add_executable(parallelrun parallelrun.cxx)
target_link_libraries(parallelrun ${LIBS_OBJ} cpuKernel)

add_executable(unpartition unpartition.cxx)
target_link_libraries(unpartition ${LIBS_OBJ} cpuKernel)

add_executable(ijparallelrun ijparallelrun.cxx)
target_link_libraries(ijparallelrun ${LIBS_OBJ} cpuKernel)

add_executable(Nparallel Nparallel.cxx)
target_link_libraries(Nparallel ${LIBS_OBJ} cpuKernel)

add_executable(skip_tree skip_tree.cxx)
target_link_libraries(skip_tree ${LIBS_OBJ} cpuKernel)

add_executable(overlap_comm overlap_comm.cxx)
target_link_libraries(overlap_comm ${LIBS_OBJ} cpuKernel)

add_definitions(-DVTK)
add_executable(vtk_test vtk.cxx ../kernel/cpuSphericalLaplace.cxx)
target_link_libraries(vtk_test ${LIBS_OBJ} )

add_executable(fdgl fdgl.cxx ../kernel/cpuSphericalLaplace.cxx ../kernel/cpuVanDerWaals.cxx ../kernel/cpuP2P.cxx)
target_link_libraries(fdgl ${LIBS_OBJ})

add_executable(poisson poisson.cxx ../kernel/cpuSphericalLaplace.cxx ../kernel/cpuVanDerWaals.cxx ../kernel/cpuP2P.cxx)
target_link_libraries(poisson ${LIBS_OBJ})
