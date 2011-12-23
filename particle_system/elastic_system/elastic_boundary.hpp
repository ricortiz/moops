#ifndef ELASTIC_BOUNDARY_HPP
#define ELASTIC_BOUNDARY_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    ElasticBoundary
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name ElasticBoundary
/// @section Description 
/// @section See also

#include <vector>
#include "geometry/surface.hpp"

template<typename spring_system_type, 
	 typename fluid_solver_type, 
	 template<typename,int,int> class integration_policy>
class ElasticBoundary : public spring_system_type, public integration_policy<ElasticBoundary<spring_system_type,fluid_solver_type,integration_policy>,5,4 >
{ 
    public:
        typedef typename spring_system_type::value_type    value_type;
        typedef typename spring_system_type::spring_type   spring_type;
        typedef typename spring_system_type::particle_type particle_type;
        typedef integration_policy<ElasticBoundary<spring_system_type,fluid_solver_type,integration_policy>,5,4 > time_integrator_type;
    private:
        size_t m_ode_size;
        fluid_solver_type                m_fluid_solver;

    public:
        ElasticBoundary(size_t ode_size) : m_ode_size(ode_size), spring_system_type(ode_size), time_integrator_type(ode_size), m_fluid_solver(ode_size/3) {  }
        ~ElasticBoundary() {}
        inline fluid_solver_type &fluid_solver() { return m_fluid_solver; }
        inline fluid_solver_type const &fluid_solver() const { return m_fluid_solver; }
        
        void run(value_type timestep)
        {
            this->integrate(timestep);
            this->time() += timestep;
        }

        void operator()(value_type time, value_type *x, value_type *v)
        {
            std::copy(x,x+m_ode_size,this->positions());
            this->update_forces(time);
            m_fluid_solver(x,v,this->forces());
        }

        void Explicit(value_type time, const value_type *x, value_type *v)
        {
            std::copy(x,x+m_ode_size,this->positions());
            this->update_forces(time);
            m_fluid_solver.ExplicitOperator(x,v,this->forces());
        }
        void Implicit(value_type time, const value_type *x, value_type *v)
        {
            std::copy(x,x+m_ode_size,this->positions());
            this->update_forces(time);
            m_fluid_solver.ImplicitOperator(x,v,this->forces());
        }
};

template< typename _spring_system_type, typename _fluid_solver_type, template<typename,int,int> class _integration_policy>
struct immersed_structure_traits<ElasticBoundary<_spring_system_type,_fluid_solver_type,_integration_policy> >
{
    typedef typename _spring_system_type::value_type    value_type;
};

#endif
