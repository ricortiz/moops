#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "particle_system/elastic_system/elastic_boundary.hpp"
#include "particle_system/particle_system.hpp"
#include "particle_system/time_integrator.hpp"
#include "particle_system/fluid_solver.hpp"

template<typename Derived>
class Surface : public ParticleSystem<Derived>,
                public ElasticBoundary<Derived>,
                public TimeIntegrator<Derived>,
                public FluidSolver<Derived>
{
    
    public:
        Surface(size_t num_particles) : ParticleSystem<Derived>(num_particles), TimeIntegrator<Derived>(3*num_particles) {}
        ~Surface() {}

        inline Derived &derived()
        {
            return *static_cast<Derived*>(this);
        }

        template<typename value_type>
        inline void run(value_type timestep)
        {
            integrate(this->time(),timestep);
            this->time() += timestep;
        }

        
};


#endif
