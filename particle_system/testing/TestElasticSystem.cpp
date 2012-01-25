#include<cassert>
#include<list>
#include<vector>
#include<map>
#include "utils/logger.hpp"

// Logger logger;

template<typename T>
struct surface_traits;

#include "math/fluid_solver/stokes/cpu_stokes_solver.hpp"
#include "math/ode_solver/euler/forward_euler.hpp"
#include "math/ode_solver/euler/backward_euler.hpp"
#include "math/ode_solver/sdc/explicit_sdc.hpp"
#include "math/ode_solver/sdc/semi_implicit_sdc.hpp"
#include "examples/valveless_heart/valveless_heart.hpp"
#include "geometry/surface.hpp"

int main()
{
    typedef CpuStokesSolver<double> fluid_solver;
    typedef ForwardEuler fetime_integrator;
    typedef BackwardEuler<double> betime_integrator;
    typedef ExplicitSDC<double> esdc_integrator;
    typedef SemiImplicitSDC<double> sisdc_integrator;
    typedef HeartPump<double,fluid_solver,esdc_integrator> heart_pump;

    heart_pump pump(20,200);
    pump.fluid_solver().setDelta(.01);
    pump.fluid_solver().initMaps(pump.elasticBoundary());

    pump.run(.01);
    pump.run(.01);

}

