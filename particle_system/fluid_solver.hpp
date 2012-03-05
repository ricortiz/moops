#ifndef FLUID_SOLVER_HPP
#define FLUID_SOLVER_HPP
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

template<typename Derived>
class FluidSolver
{
    protected:
        typedef typename Traits<Derived>::fluid_solver_type fluid_solver_type;

    private:
        fluid_solver_type m_fluid_solver;

    public:
        FluidSolver(size_t num_particles) : m_fluid_solver(num_particles) {}
        
        inline Derived &derived()
        {
            return *static_cast<Derived*>(this);
        }

        template<typename value_type>
        inline void operator()(value_type t, value_type *x, value_type *v)
        {
            derived().computeForces(t);
            m_fluid_solver(t,x,v,derived().forces());
        }

        template<typename value_type>
        inline void operator()(value_type t, value_type *x, value_type *v, size_t num_targets)
        {
            m_fluid_solver(t,x,v,derived().positions(),derived().forces(),num_targets);
        }
        
        template<typename value_type>
        inline void Explicit(value_type t, value_type *x, value_type *v)
        {
            derived().computeForces(t);
            m_fluid_solver.Explicit(t,x,v,derived().forces());
        }

        template<typename value_type>
        inline void Implicit(value_type t, value_type *x, value_type *v)
        {
            derived().computeForces(t);
            m_fluid_solver.Implicit(t,x,v,derived().forces());
        }
        fluid_solver_type &fluid_solver() { return m_fluid_solver; }
};



#endif
