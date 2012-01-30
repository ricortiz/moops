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

class ForwardEuler
{

    size_t ode_size;
    public:
        ForwardEuler(size_t _ode_size) : ode_size(_ode_size) {}
        
        template<typename function_type, typename value_type>
        inline void operator()(function_type &F, value_type t, value_type *x, value_type *xold, value_type *v, value_type dt)
        {
            F(t, xold, v);
            for(size_t i = 0; i < ode_size; ++i)
                x[i] = xold[i] + dt * v[i];
        }

        template<typename function_type, typename value_type>
        inline void operator()(function_type &, value_type, value_type *x, value_type *xold, value_type *v, value_type *vold, value_type dt)
        {
            for(size_t i = 0; i < ode_size; ++i)
                x[i] = xold[i] + dt * vold[i];
        }

        template<typename function_type, typename value_type>
        inline void operator()(function_type &F, value_type t, value_type *x, value_type *v, value_type dt)
        {
            for(size_t i = 0; i < ode_size; ++i)
                x[i] += dt * v[i];
            F(t+dt,x,v);
        }
        
        template<typename value_type>
        inline void operator()(value_type *x, const value_type *v, value_type dt)
        {
            for(size_t i = 0; i < ode_size; ++i)
                x[i] += dt * v[i];
        }

        template<typename value_type>
        inline void operator()(value_type *xnew, value_type *xold, const value_type *v, value_type dt)
        {
            for(size_t i = 0; i < ode_size; ++i)
                xnew[i] = xold[i] + dt * v[i];
        }

};



#endif
