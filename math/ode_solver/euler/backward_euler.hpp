#ifndef BACKWARD_EULER_HPP
#define BACKWARD_EULER_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    BackwardEuler
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================

#include "math/ode_solver/euler/forward_euler.hpp"
#include "math/nonlinear_solver/inexact_newton.hpp"

template<typename value_type, typename function_type, int gmres_iterations = 40, int gmres_restarts = 10>
class BackwardEuler
{
    protected:
        typedef InexactNewtonMethod<value_type, gmres_iterations, gmres_restarts> newton_solver_type;
        typedef ForwardEuler<value_type, function_type>                            forward_euler_solver_type;
        class implicit_operator
        {
                value_type m_t;
                value_type m_dt;
                std::vector<value_type> m_b;
                std::vector<value_type> m_Vx;
                function_type &m_V;

            public:

                implicit_operator ( function_type &V ) : m_V ( V ), m_b ( m_ode_size ) {}

                inline void operator() ( const value_type *x, value_type *Fx )
                {
                    std::fill ( m_Vx.begin(), m_Vx.end(), 0.0 );
                    m_V ( m_t, x, &m_Vx[0] );
                    for ( size_t i = 0; i < m_ode_size; ++i )
                        Fx[i] = x[i] - m_dt * m_Vx[i] - m_b[i];
                }

                void setParameters ( value_type t, value_type *x, const value_type *v, value_type dt, bool with_predictor = true )
                {
                    m_t = t;
                    m_dt = dt;
                    std::copy ( x, x + m_ode_size, m_b.begin() );
                    if ( with_predictor )
                        m_forward_euler_solver ( &m_b[0], v, dt );
                }

                const value_type *getRHS() const { return &m_b[0]; }
        };
        
    private:
        implicit_operator m_F;
        forward_euler_solver_type m_forward_euler_solver;
        newton_solver_type m_newton_solver;
        size_t m_ode_size;

    public:
        BackwardEuler (  function_type &F ) : m_F ( F ), m_forward_euler_solver ( F ), m_ode_size ( F.data_size() ), m_newton_solver ( F.data_size() ) {}

        inline void operator() ( value_type t, value_type *x, value_type *v, value_type dt )
        {
            m_F.setParameters ( t, x, v, dt );
            m_newton_solver ( m_F, x, 1e-6, 1e-6 );
            value_type inv_dt = 1.0 / dt;
            value_type *b = m_F.getRhs();
            for ( size_t i = 0; i < m_ode_size; ++i )
                v[i] = inv_dt * ( x[i] - b[i] );
        }

        inline void operator() ( value_type t, value_type *x, const value_type *xold, value_type *v, value_type dt )
        {
            m_F.setParameters ( t, xold, v, dt );
            m_newton_solver ( m_F, x, 1e-6, 1e-6 );
            value_type inv_dt = 1.0 / dt;
            value_type *b = m_F.getRhs();
            for ( size_t i = 0; i < m_ode_size; ++i )
                v[i] = inv_dt * ( x[i] - b[i] );
        }

        inline void operator() ( value_type t, value_type *x, const value_type *xold, value_type *v, const value_type *vold, value_type dt )
        {
            m_F.setParameters ( t, xold, vold, dt );
            m_newton_solver ( m_F, x, 1e-6, 1e-6 );
            value_type inv_dt = 1.0 / dt;
            value_type *b = m_F.getRhs();
            for ( size_t i = 0; i < m_ode_size; ++i )
                v[i] = inv_dt * ( x[i] - b[i] );
        }


};

#endif
