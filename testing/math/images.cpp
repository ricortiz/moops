#include <iostream>
#include<vector>
#include<algorithm>
#include<iterator>
#include<cstdlib>

#include "math/fluid_solver/stokes/nbody_cpu/cpu_compute_velocity.hpp"

template<typename value_type>
struct random_generator
{
    value_type operator()()
    {
        return rand()/((value_type)RAND_MAX+1);
    }
};

int images(int ,char **)
{
    srand(0);
    typedef double value_type;
    size_t num_sources = 100;
    size_t num_targets = 100;
    size_t size_sources = 3*num_sources;
    size_t size_targets = 3*num_targets;
    std::vector<value_type> sources(size_sources), targets(3), velocities(3,0), forces(size_sources);
    value_type delta = .01;
    
    std::generate(sources.begin(),sources.end(),random_generator<value_type>());
    std::generate(targets.begin(),targets.end(),random_generator<value_type>());
    std::generate(forces.begin(),forces.end(),random_generator<value_type>());
    targets[0] = 0; targets[1] = 0; targets[2] = 0;
    for(size_t i = 0; i < 1; i+=3)
        for(size_t j = 0; j < size_sources; j+=3)
        {
            computeStokeslet(&targets[i],&velocities[i],&sources[j],&forces[j],delta);
            computeImage(&targets[i],&velocities[i],&sources[j],&forces[j],delta);
        }
    std::cout.precision(16);
    std::cout << "sources = ["; std::copy(sources.begin(), sources.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "];" << std::endl;
    std::cout << "targets = ["; std::copy(targets.begin(), targets.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "];" << std::endl;
    std::cout << "forces = ["; std::copy(forces.begin(), forces.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "];" << std::endl;
    std::cout << "velocities = ["; std::copy(velocities.begin(), velocities.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "];" << std::endl;
    return 0;
}