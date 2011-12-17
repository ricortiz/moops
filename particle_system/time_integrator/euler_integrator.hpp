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
        typedef typename immersed_structure_traits<boundary_type>::value_type          value_type;

    private:
        ForwardEuler<value_type> euler_integrator;

    public:

        inline boundary_type &boundary()
        {
            return *static_cast<boundary_type*>(this);
        }

        template<typename value_type>
        inline void integrate(value_type timestep)
        {
            value_type time = boundary().time();
            euler_integrator(boundary(),time,timestep);
        }

};



#endif
