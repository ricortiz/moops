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

template<typename value_type, typename function_type>
struct BackwardEulerFunction
{
    function_type &F;
    value_type *rhs;
    value_type t;
    value_type dt;
    size_t ode_size;

    BackwardEulerFunction(function_type &_F, value_type *_rhs, value_type _t, value_type _dt, size_t _ode_size ) : F(_F), rhs(_rhs), t(_t), dt(_dt), ode_size(_ode_size) {}

    inline void operator()(value_type *x, value_type *Fx)
    {
        std::fill(Fx, Fx + ode_size, 0.0);
        F(t, x, Fx);
        for (size_t i = 0; i < ode_size; ++i)
        {
            Fx[i] *= -dt;
            Fx[i] += x[i] - rhs[i];
        }
    }
};

template < typename value_type, int gmres_iterations = 1000, int gmres_restarts = 100 >
class BackwardEuler
{
    protected:
        typedef InexactNewtonMethod<value_type, gmres_iterations, gmres_restarts> newton_solver_type;
        typedef ForwardEuler                            forward_euler_solver_type;

    protected:
        forward_euler_solver_type forward_euler;
        newton_solver_type newton_solver;

    private:
        value_type m_t;
        value_type m_dt;

        size_t ode_size;

    public:
        BackwardEuler(size_t _ode_size) : forward_euler(_ode_size), ode_size(_ode_size), newton_solver(_ode_size)  {}
        
        template<typename function_type>
        inline void operator()(function_type F, value_type t, value_type *x, value_type *v, value_type dt)
        {
	  operator()(F,t,x,x,v,v,dt);
        }

        template<typename function_type>
        inline void operator()(function_type F, value_type t, value_type *x, value_type *xold, value_type *v, value_type dt)
        {
            BackwardEulerFunction<value_type, function_type> G(F, xold, t, dt, ode_size);
            newton_solver(G, x, 1e-15, 1e-15);
            value_type inv_dt = 1.0 / dt;
            for (size_t i = 0; i < ode_size; ++i)
                v[i] = inv_dt * (x[i] - xold[i])-newton_solver.f()[i];
        }

};

#endif
