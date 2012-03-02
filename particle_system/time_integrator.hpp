#ifndef TIME_INTEGRATOR_HPP
#define TIME_INTEGRATOR_HPP
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
class TimeIntegrator
{
    protected:
        typedef typename Traits<Derived>::time_integrator_type time_integrator_type;

    private:
        time_integrator_type time_integrator;

    public:
        TimeIntegrator(size_t ode_size) : time_integrator(ode_size) {}
        
        inline Derived &derived()
        {
            return *static_cast<Derived*>(this);
        }

        template<typename value_type>
        inline void integrate(value_type t, value_type timestep)
        {
            time_integrator(derived(),t,derived().positions(),derived().velocities(),timestep);
        }
        
};



#endif
