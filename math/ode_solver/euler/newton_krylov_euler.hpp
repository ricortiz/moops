#ifndef NEWTON_KRYLOV_EULER_HPP
#define NEWTON_KRYLOV_EULER_HPP
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

#include "math/ode_solver/euler/backward_euler.hpp"
#include "math/nonlinear_solver/inexact_newton.hpp"

template <typename value_type, typename ode_rhs_type, int gmres_iterations = 100, int gmres_restarts = 10>
class NewtonKrylovEuler
{
    protected:
        typedef InexactNewtonMethod<value_type, gmres_iterations, gmres_restarts> newton_solver_type;
        class nonlinear_operator
        {
            protected:
                typedef BackwardEuler<value_type, ode_rhs_type> euler_type;

            private:
                euler_type m_euler;
                value_type m_t;
                value_type m_dt;
                ode_rhs_type &m_V;
                value_type *m_v;
                std::vector<value_type> m_x;

            public:

                nonlinear_operator(ode_rhs_type &V) : m_euler(V), m_V(V), m_x(V.ode_size(), 0.0)
                {

                }

                inline void operator()(value_type *x, value_type *Fx)
                {
                    std::copy(x, x + m_x.size(), m_x.begin());
                    m_V(m_t, x, m_v);
                    m_euler(m_t, x, m_v, m_dt);
                    std::transform(m_x.begin(), m_x.end(), x, Fx, std::minus<value_type>());
                }

                void init(value_type t, value_type dt, value_type *v)
                {
                    m_t = t;
                    m_dt = dt;
                    m_v = v;
                }
        };

    public:
        nonlinear_operator m_F;
        newton_solver_type m_newton_solver;

    public:
        NewtonKrylovEuler(ode_rhs_type &F) : m_F(F), m_newton_solver(F.ode_size()) {}

        inline void operator()(value_type t, value_type *x, value_type *v, value_type dt)
        {
            m_F.init(t, dt, v);
            m_newton_solver(m_F, x, 1e-12, 1e-6);
        }

};

#endif
