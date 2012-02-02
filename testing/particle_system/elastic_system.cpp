#include<cassert>
#include<list>
#include<vector>
#include<map>
#include<fstream>
#include "utils/logger.hpp"

// Logger logger;

template<typename T>
struct Traits;

#include "math/fluid_solver/stokes/cpu_stokes_solver.hpp"
#include "math/ode_solver/euler/forward_euler.hpp"
#include "math/ode_solver/euler/backward_euler.hpp"
#include "math/ode_solver/sdc/explicit_sdc.hpp"
#include "math/ode_solver/sdc/semi_implicit_sdc.hpp"
// #include "examples/valveless_heart/valveless_heart.hpp"
#include "examples/swarm/swarm.hpp"

int elastic_system(int, char **)
{
  
    typedef CpuStokesSolver<double> fluid_solver;
    typedef ForwardEuler fetime_integrator;
    typedef BackwardEuler<double> betime_integrator;
    typedef ExplicitSDC<double> esdc_integrator;
    typedef SemiImplicitSDC<double> sisdc_integrator;
//     typedef HeartPump<double,fluid_solver,esdc_integrator> heart_pump;
    typedef Swarm<double,fluid_solver,fetime_integrator> swarm_type;

    swarm_type swarm(6,100,12,21,12100);
    std::ofstream output("swimmers_10M.dat");
    output.precision(16);
    swarm.print_positions(output,true);
    output.flush();
    output.close();
//     swarm.print_springs(std::cout);
//     swarm.print_tail_springs(std::cout);
    return 0;
}

