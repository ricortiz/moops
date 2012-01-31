#include <iostream>
#include<vector>
#include<algorithm>
#include<cstdlib>

#include "math/fluid_solver/stokes/cpu_stokes_solver.hpp"

template<typename value_type>
struct random_generator
{
    value_type operator()()
    {
        return rand()/((value_type)RAND_MAX+1);
    }
};

int cpu_stokes_solver(int ,char **)
{
    srand(0);
    typedef double value_type;
    size_t num_sources = 100;
    size_t num_targets = 100;
    size_t size_sources = 3*num_sources;
    size_t size_targets = 3*num_targets;
    std::vector<value_type> sources(size_sources), targets(size_targets), velocities(size_targets), forces(size_sources);
    value_type delta = .01;
    CpuStokesSolver<double> solver(0);
    solver.setDelta(delta);

    std::generate(sources.begin(),sources.end(),random_generator<value_type>());
    std::generate(targets.begin(),targets.end(),random_generator<value_type>());
    std::generate(velocities.begin(),velocities.end(),random_generator<value_type>());
    std::generate(forces.begin(),forces.end(),random_generator<value_type>());

    solver(0,&targets[0],&velocities[0],&sources[0],&forces[0],num_sources,num_targets);
    solver(0,&targets[0],&velocities[0],&targets[0],&forces[0],num_sources,num_targets);
    solver(0,&sources[0],&velocities[0],&sources[0],&forces[0],num_sources,num_targets);
    
    return 0;
}