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

template<typename value_type, typename function_type>
class ForwardEuler
{
    function_type &m_F;
    size_t m_ode_size;
    
    public:
        ForwardEuler(function_type &F) : m_F(F), m_ode_size(F.ode_size()) {}
        
        inline void operator()(value_type t, value_type *x, value_type *v, value_type dt)
        {
            m_F(t,x,v);
            for (size_t i = 0; i < m_ode_size; ++i)
                x[i] += dt*v[i];
        }
        
        inline void operator()(value_type *x, const value_type *v, value_type dt)
        {
            for (size_t i = 0; i < m_ode_size; ++i)
                x[i] += dt*v[i];
        }
        
        inline void operator()(value_type *x, const value_type *xold, const value_type *v, value_type dt)
        {
            for (size_t i = 0; i < m_ode_size; ++i)
                x[i] = xold[i] + dt*v[i];
        }

};



#endif
