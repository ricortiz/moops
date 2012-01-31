#ifndef TYPE_BINDER_HPP
#define TYPE_BINDER_HPP

#include<list>
#include<cmath>
#include<cassert>
#include<sstream>
#include<algorithm>
#include<iterator>

#ifdef USE_QT_GUI
#include <QVTKApplication.h>
#include "gui/gui.hpp"
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
#endif

template<typename _value_type>
struct TypeBinder
{
    typedef TypeBinder<_value_type>                             Types;
    typedef _value_type value_type;
    typedef Particle<value_type>                                particle_type;
#ifdef USE_CUDA_FLUID_SOLVER
    typedef GpuStokesSolver<float>                              fluid_solver;
#else
    typedef CpuStokesSolver<value_type>                         fluid_solver;
#endif

    typedef ForwardEuler                                        forward_euler;
    typedef BackwardEuler<value_type>                           backward_euler;

    enum
    {
        sdc_nodes = 3,
        sdc_corrections = 3
    };
    typedef Integrator<value_type, gauss_lobatto, sdc_nodes> spectral_integrator;

    typedef ExplicitSDC<value_type, spectral_integrator, sdc_corrections>                             explicit_sdc;
    typedef SemiImplicitSDC<value_type, spectral_integrator, sdc_corrections>                         implicit_sdc;

    typedef implicit_sdc                                        time_integrator;

    typedef HeartPump<value_type, fluid_solver, time_integrator>  heart_pump_surface;
    typedef Swarm<value_type, fluid_solver, time_integrator>      swarm_surface;

    typedef vtkStorageWrapper<swarm_surface>                            vtk_storage;
    typedef IO::VtkWriter<vtk_storage, vtkXMLPolyDataWriter>                 vtk_writer;


};

#endif
