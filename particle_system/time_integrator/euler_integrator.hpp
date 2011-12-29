#ifndef EULER_INTEGRATOR_HPP
#define EULER_INTEGRATOR_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    EulerIntegrator
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name SDCIntegrator -- Wrapper for the and Euler-based time integrator
/// @section Description EulerIntegrator wraps the integrator so it can be used by the simulator
///                      It contains the main method that drives the simulation: void integrate()
/// @section See also ForwardEuler

#include "math/ode_solver/euler/forward_euler.hpp"

template<typename boundary_type>
class EulerIntegrator
{
    protected:
        typedef EulerIntegrator<boundary_type>                                         self_type;
        typedef typename immersed_structure_traits<boundary_type>::fluid_solver_type       fluid_solver_type;
        typedef typename immersed_structure_traits<boundary_type>::value_type          value_type;

    private:
        size_t    m_ode_size;
        ForwardEuler<value_type,boundary_type> m_euler_integrator;
	fluid_solver_type &m_fluid_solver;

    public:
        
        EulerIntegrator(size_t ode_size) : m_ode_size(ode_size), m_euler_integrator(derived()), m_fluid_solver(derived().fluid_solver()) {}
        
        inline boundary_type &derived()
        {
            return *static_cast<boundary_type*>(this);
        }

        template<typename value_type>
        inline void integrate(value_type timestep)
        {
            value_type time = derived().time();
            value_type positions = derived().positions();
            value_type velocities = derived().velocities();
            size_t data_size = derived().data_size();
            m_euler_integrator(time,positions,velocities,timestep,data_size);
        }       
        
        void operator()(value_type time, value_type *x, value_type *v)
        {
            std::copy(x, x + m_ode_size, derived().positions());
            derived().update_forces(time);
            m_fluid_solver(x, v, derived().forces());
        }

        size_t ode_size() { return derived().data_size(); }
        
};



#endif
