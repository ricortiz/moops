#include <iostream>
#include<vector>
#include<algorithm>
#include<iterator>
#include<cstdlib>

#include "utils/logger.hpp"

Logger logger;
#include "math/fluid_solver/stokes/hybrid_fmm_stokes_solver.hpp"

template<typename value_type>
struct random_generator
{
    value_type operator()()
    {
        return rand()/((value_type)RAND_MAX+1);
    }
};

int hybrid_fmm_stokes_solver(int ,char **)
{
    srand(0);
    typedef float value_type;
    size_t num_sources = 1 << 7;
    size_t num_targets = 1 << 7;
    size_t size_sources = 3*num_sources;
    size_t size_targets = 3*num_targets;
    std::vector<value_type> sources(size_sources), targets(size_targets), velocities(size_targets), forces(size_sources);
    value_type delta = .01;
    value_type domain[2][3] = {{0,0,0},{10,10,10}};
    HybridFmmStokesSolver<value_type> fmm(num_sources);
    std::generate(sources.begin(),sources.end(),random_generator<value_type>());
    std::generate(targets.begin(),targets.end(),random_generator<value_type>());
    std::generate(forces.begin(),forces.end(),random_generator<value_type>());    
    std::copy(sources.begin(),sources.end(),std::ostream_iterator<value_type>(std::cout," ")); std::cout << std::endl;    
    std::copy(targets.begin(),targets.end(),std::ostream_iterator<value_type>(std::cout," ")); std::cout << std::endl;    
    std::copy(forces.begin(),forces.end(),std::ostream_iterator<value_type>(std::cout," ")); std::cout << std::endl;
    fmm.setDelta(delta);
    fmm(0,&sources[0],&velocities[0],&forces[0]);
    logger.printTimer();
    return 0;
}