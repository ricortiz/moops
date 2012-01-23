#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "particle_system/elastic_system/elastic_boundary.hpp"
#include "particle_system/particle_system.hpp"
#include "particle_system/time_integrator.hpp"

template<typename T> struct surface_traits;

template<typename Derived>
class Surface : ElasticBoundary<Derived>, ParticleSystem<Derived>, TimeIntegrator<Derived>
{
    protected:
        typedef typename surface_traits<Derived>::fluid_solver_type fluid_solver_type;
        typedef typename surface_traits<Derived>::time_integrator time_integrator_type;
        typedef typename surface_traits<Derived>::value_type value_type;

    protected:
	fluid_solver_type m_fluid_solver;
	time_integrator_type m_time_integrator;
	
    public:
        Surface() {}
        ~Surface() {}

        inline Derived *derived()
        {
            return static_cast<Derived*>(this);
        }

        inline void eval(value_type time, value_type time_step)
	{
	  m_time_integrator();
	}
	inline void eval(value_type t, value_type *x, value_type *v)
	{
            derived().updateForces(time);
            m_fluid_solver(x, v, this->forces());    
	}
        inline void run(value_type timestep)
        {
            this->integrate(this->time(), timestep);
            this->time() += timestep;
        }
};


#endif
