#ifndef FORWARD_EULER_HPP
#define FORWARD_EULER_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    ForwardEuler
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================

template<typename value_type>
class ForwardEuler
{
    public:
        template<typename function_type>
        inline void operator()(function_type &V, value_type t, const value_type dt)
        {
            value_type *x = V.positions();
            value_type *v = V.velocities();

            size_t ode_size = 3*V.particles_size();
            V(t,x,v);
            for (size_t i = 0; i < ode_size; ++i)
                x[i] += dt*v[i];
        }

};



#endif
