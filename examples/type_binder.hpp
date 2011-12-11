#ifndef TYPE_BINDER_HPP
#define TYPE_BINDER_HPP

#ifdef USE_QT_GUI
#include <QVTKApplication.h>
#include "gui/gui.hpp"
#endif

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkXMLPolyDataWriter.h>

#include "particle_system/time_integrator/sdc_integrator.hpp"
#include "particle_system/time_integrator/euler_integrator.hpp"
#include "particle_system/elastic_system/elastic_boundary.hpp"
#include "particle_system/elastic_system/spring_system.hpp"
#include "particle_system/particle_system.hpp"
#include "particle_system/fluid_solver/fmm_stokes_solver.hpp"
#ifdef USE_CUDA_FLUID_SOLVER
#include "particle_system/fluid_solver/cuda_stokes_solver.hpp"
#endif
#include "particle_system/fluid_solver/direct_stokes_solver.hpp"
#include "math/ode_solver/sdc/explicit_sdc.hpp"
#include "math/ode_solver/euler/forward_euler.hpp"
#include "geometry/torus_geometry.hpp"
#include "geometry/oval_geometry.hpp"

#include "io/write_vtu.hpp"
#include "particle_system/particle_markers.hpp"
#include "geometry/sine_geometry.hpp"
#include "swarm/swarm.hpp"
#include "valveless_heart/valveless_heart.hpp"

template<typename value_type, int num_particles, int sdc_nodes = 5, int data_size= 3*num_particles, int fmm_max_particles = 0, int fmm_order = 0>
struct TypeBinder
{
    typedef TypeBinder<value_type,num_particles,sdc_nodes,data_size,fmm_max_particles,fmm_order> Types;
    
    typedef ClenshawCurtis<value_type,data_size,sdc_nodes>                      sdc_integrator_type;
    typedef ExplicitSDC<value_type,sdc_integrator_type,sdc_nodes,sdc_nodes>     sdc_type;
    typedef ForwardEuler<value_type>                                            euler_type;
    typedef Particle<value_type>                                                particle_type;
#ifdef USE_CUDA_FLUID_SOLVER
    typedef vtkParticleSystemStorage<value_type,particle_type,PSYS::SURFACE,vtkFloatArray> surface_storage_type;
    typedef vtkParticleSystemStorage<value_type,particle_type,PSYS::VOLUME,vtkFloatArray> volume_storage_type;
    typedef ParticleSystem<value_type,sdc_type,num_particles,PSYS::SURFACE,surface_storage_type> particle_system_type;
    typedef ParticleSystem<value_type,euler_type,num_particles,PSYS::VOLUME,volume_storage_type> particle_system_tracers_type;
#else
    typedef ParticleSystem<value_type,sdc_type,num_particles>                   particle_system_type;
    typedef ParticleSystem<value_type,euler_type,num_particles,PSYS::VOLUME>    particle_system_tracers_type;
#endif
    typedef HeartPump<value_type>                                               heart_pump_surface_type;
    typedef Swarm<value_type>                                                   swarm_surface_type;                                              
    typedef SineGeometry<value_type>                                            sine_type;
    typedef OvalGeometry<value_type>                                            oval_type;
    typedef SpringSystem<heart_pump_surface_type,particle_system_type>          heart_pump_spring_system_type;
    typedef SpringSystem<swarm_surface_type,particle_system_type>               swarm_spring_system_type;
#ifdef USE_CUDA_FLUID_SOLVER
    typedef CudaStokesSolver<particle_system_type>                              fluid_solver_type;
#else
    typedef DirectStokesSolver<particle_system_type>                            fluid_solver_type;
    typedef FMMStokesSolver<particle_system_type,fmm_max_particles,fmm_order> fmm_fluid_solver_type;
#endif
    typedef ElasticBoundary<heart_pump_spring_system_type,fluid_solver_type,SDCIntegrator> heart_pump_boundary_type;
    typedef ElasticBoundary<swarm_spring_system_type,fluid_solver_type,SDCIntegrator> swarm_boundary_type;
    typedef ParticleMarkers<particle_system_tracers_type,heart_pump_boundary_type,fluid_solver_type,EulerIntegrator> heart_pump_volume_type;
    typedef ParticleMarkers<particle_system_tracers_type,swarm_boundary_type,fluid_solver_type,EulerIntegrator> swarm_volume_type;
    
    typedef IO::VTKWriter<vtkXMLPolyDataWriter>                                 vtk_writer;
};

#endif
