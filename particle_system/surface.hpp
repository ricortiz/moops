#ifndef SURFACE_HPP
#define SURFACE_HPP
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
#include "particle_system/elastic_system/elastic_boundary.hpp"
#include "particle_system/particle_system.hpp"
#include "particle_system/time_integrator.hpp"
#include "particle_system/fluid_solver.hpp"

template<typename Derived>
class Surface : public ParticleSystem<Derived>,
                public ElasticBoundary<Derived>,
                public TimeIntegrator<Derived>,
                public FluidSolver<Derived>
{
    
    public:
        Surface(size_t num_particles) : ParticleSystem<Derived>(num_particles), TimeIntegrator<Derived>(3*num_particles), FluidSolver<Derived>(num_particles) {}
        ~Surface() {}

        inline Derived &derived()
        {
            return *static_cast<Derived*>(this);
        }

        inline ParticleSystem<Derived> &particleSystem()
        {
            return *static_cast<ParticleSystem<Derived>*>(this);
        }
        
        inline ElasticBoundary<Derived> &elasticBoundary()
        {
            return *static_cast<ElasticBoundary<Derived>*>(this);
        }

        inline SpringSystem<Derived> &springSystem()
        {
            return *static_cast<SpringSystem<Derived>*>(this);
        }
        
        inline TimeIntegrator<Derived> &timeIntegrator()
        {
            return *static_cast<TimeIntegrator<Derived>*>(this);
        }
        
        inline FluidSolver<Derived> &fluidSolver()
        {
            return *static_cast<FluidSolver<Derived>*>(this);
        }
        
        template<typename value_type>
        inline void run(value_type timestep)
        {
            this->integrate(this->time(),timestep);
            this->time() += timestep;
        }

        
};


#endif
