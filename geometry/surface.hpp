#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "particle_system/elastic_system/elastic_boundary.hpp"
#include "particle_system/particle_system.hpp"

template<typename T> struct surface_traits;

template<typename Derived>
class Surface : ElasticBoundary<Derived>, ParticleSystem<Derived>
{
    protected:
        typedef typename surface_traits<Derived>::fluid_solver_type fluid_solver_type;
        typedef typename surface_traits<Derived>::value_type value_type;

    protected:
	fluid_solver_type m_fluid_solver;
	
    public:
        Surface() {}
        ~Surface() {}

        inline Derived *derived()
        {
            return static_cast<Derived*>(this);
        }

	void eval(value_type t, value_type *x, value_type *v)
	{
            logger.startTimer("odeRHSeval");
            derived().updateForces(time);
            m_fluid_solver(x, v, derived().forces());
            logger.stopTimer("odeRHSeval");        
	}

};


#endif
