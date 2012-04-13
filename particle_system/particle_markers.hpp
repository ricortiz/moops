#ifndef PARTICLE_MARKERS_HPP
#define PARTICLE_MARKERS_HPP
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

#include "particle_system.hpp"
#include "storage/particle_system_storage.hpp"

/** \class ParticleMarkers
 *  \ingroup ParticleSystem_Module
 *  \brief This is a helper class that keep tracks of markers in the particle system
 *  \tparam surface_type is the surface particle system associated to this markers,  it defines the fluid solver used here.
 *  \tparam time_integrator is the time stepper used to move the markers, usually inexpensive and not necesarily the same used by the surface_type
 */
template<typename surface_type, typename time_integrator>
class ParticleMarkers : public ParticleSystem<ParticleMarkers<surface_type, time_integrator> >
{
    protected:
        typedef ParticleSystem<ParticleMarkers<surface_type, time_integrator> > base_type;

    protected:
        time_integrator m_integrator;

    public:
        ParticleMarkers(size_t num_particles) : base_type(num_particles), m_integrator(3*num_particles) {  }

        template<typename value_type>
        void run(value_type timestep, surface_type &surface)
        {
            size_t num_targets = this->particles_size();
            surface(this->time(), this->positions(), this->velocities(), this->particles_size());
            m_integrator(this->positions(), this->velocities(), timestep);
            this->time() += timestep;
        }
};


template<typename surface_type, typename time_integrator>
struct Traits<ParticleMarkers<surface_type, time_integrator> >
{
    typedef typename Traits<surface_type>::value_type value_type;
    typedef typename Traits<surface_type>::particle_type particle_type;
    typedef ParticleSystemStorage< value_type, particle_type, VOLUME > storage_type;
};

#endif
