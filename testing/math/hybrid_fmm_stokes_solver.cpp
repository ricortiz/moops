#include <iostream>
#include<vector>
#include<algorithm>
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
    typedef double value_type;
    size_t num_sources = 100;
    size_t num_targets = 100;
    size_t size_sources = 3*num_sources;
    size_t size_targets = 3*num_targets;
    std::vector<value_type> sources(size_sources), targets(size_targets), velocities(size_targets), forces(size_sources);
    value_type delta = .01;
    value_type domain[2][3] = {{0,0,0},{10,10,10}};
    HybridFmmStokesSolver<value_type> fmm(num_sources);
    fmm.setDomain(domain);
    fmm.setDelta(delta);
    
    return 0;
}