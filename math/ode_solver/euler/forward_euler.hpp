#ifndef FORWARD_EULER_HPP
#define FORWARD_EULER_HPP
/****************************************************************************
** MOOPS -- Modular Object Oriented Particle Simulator
** Copyright (C) 2011-2012  Ricardo Ortiz <ortiz@unc.edu>
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/

class ForwardEuler
{

    size_t ode_size;
    public:
        ForwardEuler(size_t _ode_size) : ode_size(_ode_size) {}
        
        template<typename function_type, typename value_type>
        inline void operator()(function_type &F, value_type t, value_type *x, value_type *xold, value_type *v, value_type dt)
        {
            F(t, xold, v);
            size_t i;
            #pragma omp parallel for private(i)
            for(i = 0; i < ode_size; ++i)
                x[i] = xold[i] + dt * v[i];
        }

        template<typename function_type, typename value_type>
        inline void operator()(function_type &, value_type, value_type *x, value_type *xold, value_type *v, value_type *vold, value_type dt)
        {
            size_t i;
            #pragma omp parallel for private(i)
            for(i = 0; i < ode_size; ++i)
                x[i] = xold[i] + dt * vold[i];
        }

        template<typename function_type, typename value_type>
        inline void operator()(function_type &F, value_type t, value_type *x, value_type *v, value_type dt)
        {
            size_t i;
            #pragma omp parallel for private(i)
            for(i = 0; i < ode_size; ++i)
                x[i] += dt * v[i];
            F(t+dt,x,v);
        }
        
        template<typename value_type>
        inline void operator()(value_type *x, const value_type *v, value_type dt)
        {
            size_t i;
            #pragma omp parallel for private(i)
            for(i = 0; i < ode_size; ++i)
                x[i] += dt * v[i];
        }

        template<typename value_type>
        inline void operator()(value_type *xnew, value_type *xold, const value_type *v, value_type dt)
        {
            size_t i;
            #pragma omp parallel for private(i)
            for(size_t i = 0; i < ode_size; ++i)
                xnew[i] = xold[i] + dt * v[i];
        }

};



#endif
