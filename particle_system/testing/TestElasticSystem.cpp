#include<cassert>
#include "utils/logger.hpp"

Logger logger;

#include "particle_system/particle_system.hpp"
#include "math/fluid_solver/stokes/cpu_stokes_solver.hpp"
#include "examples/valveless_heart/valveless_heart.hpp"
#include "math/ode_solver/sdc_time_integrator.hpp"

int main()
{
  typedef ParticleSystem<double> particle_system_type;
  typedef CpuStokesSolver<double> fluid_solver_type;
  typedef SdcTimeIntegrator<double,fluid_solver_type> time_integrator_type;
  
  typedef HeartPump<particle_system_type,fluid_solver_type,time_integrator_type> heart_pump_type;
  OvalGeometry<double> oval_geometry;
  oval_geometry.setDimensions(20,200);
  oval_geometry.setX0(0,0,0);
  oval_geometry.setForcingRange(10,50);
  oval_geometry.setWaveSpeed(.001);
  heart_pump_type heart_pump(oval_geometry);
  heart_pump.run(.01);
    

}

