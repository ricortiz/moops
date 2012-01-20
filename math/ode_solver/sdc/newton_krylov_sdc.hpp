#ifndef NEWTON_KRYLOV_SDC_HPP
#define NEWTON_KRYLOV_SDC_HPP
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

#include "math/ode_solver/sdc/integrator/clenshaw_curtis.hpp"
#include "math/ode_solver/sdc/explicit_sdc.hpp"
#include "math/nonlinear_solver/inexact_newton.hpp"

template < typename value_type, typename ode_rhs_type, int gmres_iterations = 1000, int gmres_restarts = 100, int sdc_nodes = 5, int num_corrections = 5 >
class NewtonKrylovSdc
{
    protected:
        typedef InexactNewtonMethod<value_type, gmres_iterations, gmres_restarts> newton_solver_type;
        class nonlinear_operator
        {
        protected:
            typedef SDCSpectralIntegrator<value_type, 0, sdc_nodes, 2, sdc_nodes>                  spectral_integrator_type;
            typedef ExplicitSDC<value_type, ode_rhs_type, spectral_integrator_type, sdc_nodes, num_corrections>                         sdc_type;

            private:
                sdc_type m_sdc;
                value_type m_t;
                value_type m_dt;
                ode_rhs_type &m_V;

            public:

                nonlinear_operator(ode_rhs_type &V) : m_sdc(V), m_V(V)
                {
                }

                inline void operator()(value_type *x, value_type *Fx)
                {
                    m_sdc.setX0(x);
                    m_V(m_t, x, m_sdc.F(0));
                    m_sdc.predictor(m_t, m_dt);
                    m_sdc.corrector(m_t, m_dt);
                    std::transform(m_sdc.X(0),m_sdc.X(0)+m_V.ode_size(),m_sdc.X(4),Fx,std::minus<value_type>());
                }

                void init(value_type t, value_type dt, value_type *v)
                {
                    m_t = t;
                    m_dt = dt;
                    m_sdc.setF0(v);
                }
        };

    public:
        nonlinear_operator m_F;
        newton_solver_type m_newton_solver;

    public:
        NewtonKrylovSdc(ode_rhs_type &F) : m_F(F), m_newton_solver(F.ode_size()) {}

        inline void operator()(value_type t, value_type *x, value_type *v, value_type dt)
        {
            m_F.init(t, dt, v);
            m_newton_solver(m_F, x, 1e-12, 1e-6);
        }

};

#endif
