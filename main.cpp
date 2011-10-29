#include <iostream>
#include <iterator>
#include "particle_system/time_stepping/sdc_integrator.hpp"
#include "particle_system/spring_system/immersed_boundary.hpp"
#include "particle_system/spring_system/spring_system.hpp"
#include "geometry/heart_pump.hpp"
#include "geometry/swimmers.hpp"
#include "geometry/sine_geometry.hpp"
#include "particle_system/particle_system.hpp"
#include "math/ode_solver/sdc/explicit_sdc.hpp"
#include "particle_system/fluid_solver/fmm_stokes_solver.hpp"
#include "particle_system/fluid_solver/cuda_stokes_solver.hpp"
#include "geometry/torus_geometry.hpp"

// template<typename derived>
// struct class1 {};
// 
// template<template<typename> class T>
// struct class2 : public T<class2<T> >
// {
//     
// };
// 
// template<typename T>
// struct class3 {};
const size_t M = 10;
const size_t N = 100;
const size_t num_geometries = 1;
const size_t num_particles = M*N*num_geometries;
const size_t data_size = num_particles*3;
const size_t sdc_nodes = 5;
const size_t fmm_max_particles = 50;
const size_t fmm_order = 3;
typedef float value_type;

// class ode_solver :      public ExplicitSDC<value_type,ClenshawCurtis<value_type,data_size,sdc_nodes>,sdc_nodes,sdc_nodes> {};
// class particle_system : public ParticleSystem<value_type,ode_solver,num_particles> {};
// class fluid_solver :    public CudaStokesSolver<particle_system> {};
// class surface :         public HeartPump<value_type> {};
// class elastic_model :   public SpringSystem<surface,particle_system> {};
// class explicit_time_stepper_particle_simulator : public ImmersedBoundary<elastic_model,fluid_solver,SDCIntegrator> {};

typedef ClenshawCurtis<value_type,data_size,sdc_nodes> integrator_type;
typedef ExplicitSDC<value_type,integrator_type,sdc_nodes,sdc_nodes> sdc_type;
typedef ParticleSystem<value_type,sdc_type,num_particles> particle_system_type;
typedef HeartPump<value_type> heart_pump_surface_type;
typedef Swimmers<value_type> swimmers_surface_type;
typedef TorusGeometry<value_type> torus_type;
typedef SineGeometry<value_type> sine_type;
typedef SpringSystem<heart_pump_surface_type,particle_system_type> spring_system_type;
typedef FMMStokesSolver<particle_system_type,fmm_max_particles,fmm_order> fluid_solver_type;
// typedef CudaStokesSolver<particle_system_type> fluid_solver_type;
typedef ImmersedBoundary<spring_system_type,fluid_solver_type,SDCIntegrator> boundary_type;


int main(int argc, char **argv) 
{
    value_type x0[3] = {0};
    torus_type torus(x0,M,N);
    heart_pump_surface_type heart_pump;
    boundary_type boundary(.01);
    heart_pump.init(torus,boundary.positions(),boundary.particles_size());
    boundary.init(heart_pump);
    boundary.run(.01);
    value_type *p = boundary.positions();
    std::cout << "p = [";
    std::copy(p,p+data_size,std::ostream_iterator<value_type>(std::cout," "));
    std::cout << "];" <<  std::endl;
    value_type *f = boundary.forces();
    std::cout << "f = [";
    std::copy(f,f+data_size,std::ostream_iterator<value_type>(std::cout," "));
    std::cout << "];" <<  std::endl;
    value_type *v = boundary.velocities();
    std::cout << "v = [";
    std::copy(v,v+data_size,std::ostream_iterator<value_type>(std::cout," "));
    std::cout << "];" <<  std::endl;
//     std::cout << "Hello, world!" << std::endl;
    return 0;
}
