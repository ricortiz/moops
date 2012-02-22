#ifndef TYPE_BINDER_HPP
#define TYPE_BINDER_HPP

#include<list>
#include<cmath>
#include<cassert>
#include<sstream>
#include<algorithm>
#include<iterator>
#include<iterator>
#ifdef USE_QT_GUI
#include "gui/gui.hpp"
#endif
#ifdef USE_PV_COPROCESSOR
#include "utils/paraview_coprocessor.hpp"
#endif
#include "utils/logger.hpp"
#include "utils/vtk_storage_wrapper.hpp"
#include "io/write_vtu.hpp"

Logger logger;

template<typename T>
struct Traits;

#include "math/fluid_solver/stokes/cpu_stokes_solver.hpp"
#include "math/ode_solver/euler/forward_euler.hpp"
#include "math/ode_solver/euler/backward_euler.hpp"
#include "math/ode_solver/sdc/explicit_sdc.hpp"
#include "math/ode_solver/sdc/semi_implicit_sdc.hpp"
#include "examples/valveless_heart/valveless_heart.hpp"
#include "examples/swarm/swarm.hpp"

#ifdef USE_CUDA_FLUID_SOLVER
#include "math/fluid_solver/stokes/gpu_stokes_solver.hpp"
#include "math/fluid_solver/stokes/hybrid_fmm_stokes_solver.hpp"
#endif

template < typename _value_type, int _sdc_nodes = 3, int _sdc_corrections = 2 >
struct TypeBinder
{
    enum
    {
        sdc_nodes = _sdc_nodes,
        sdc_corrections = _sdc_corrections
    };
    typedef TypeBinder<_value_type>                             Types;
    typedef _value_type                                         value_type;
    typedef ParticleWrapper<value_type>                         particle_type;
#ifdef USE_CUDA_FLUID_SOLVER
    typedef GpuStokesSolver<float>                              cuda_fluid_solver;
    typedef HybridFmmStokesSolver<value_type>                   fmm_fluid_solver;
#endif
    typedef CpuStokesSolver<value_type>                         fluid_solver;

    typedef ForwardEuler                                        forward_euler;
    typedef BackwardEuler<value_type>                           backward_euler;
    typedef Integrator<value_type, gauss_lobatto, sdc_nodes>                    spectral_integrator;
    typedef ExplicitSDC<value_type, spectral_integrator, sdc_corrections>       explicit_sdc;
    typedef SemiImplicitSDC<value_type, spectral_integrator, sdc_corrections>   implicit_sdc;
    typedef explicit_sdc                                                time_integrator;
    typedef HeartPump<value_type, fluid_solver, time_integrator>        heart_pump_surface;
    typedef Swarm<value_type, fmm_fluid_solver, time_integrator>        swarm_surface;
    typedef vtkStorageWrapper<swarm_surface, vtkFloatArray>              swarm_vtk_storage;
    typedef vtkStorageWrapper<heart_pump_surface>                       heart_vtk_storage;
    typedef IO::VtkWriter<vtkXMLPolyDataWriter>                         vtk_poly_writer;
    typedef IO::VtkWriter<vtkXMLUnstructuredGridWriter>                 vtk_unstructured_writer;


};

#endif
