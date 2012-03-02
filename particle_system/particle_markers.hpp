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


template < typename particle_system_type,
typename immersed_surface_type,
template<typename> class integration_policy >
class ParticleMarkers : public particle_system_type, public integration_policy<ParticleMarkers<particle_system_type, immersed_surface_type, integration_policy> >
{
    public:
        typedef typename immersed_structure_traits<immersed_surface_type>::ode_rhs_type ode_rhs_type;
        typedef typename particle_system_type::value_type    value_type;
        typedef typename particle_system_type::particle_type particle_type;
        typedef integration_policy<ParticleMarkers<particle_system_type,immersed_surface_type,integration_policy> > time_integrator_type;

    private:
        immersed_surface_type            &m_surface;

    public:

        ParticleMarkers(immersed_surface_type &surface, size_t ode_size) : particle_system_type(ode_size), time_integrator_type(ode_size), m_surface(surface) {  }
        ~ParticleMarkers() {}

        void run(value_type timestep)
        {
            this->integrate(timestep);
            this->time() += timestep;
        }

        void operator()(value_type time, value_type *x, value_type *v)
        {
            value_type *f = m_surface.forces();
            value_type *y = m_surface.positions();
            size_t num_targets = this->particles_size();
            size_t num_sources = m_surface.particles_size();
            this->ode_rhs()(x, v, y, f, num_sources, num_targets);
        }
};

template< typename _particle_system_type, typename _immersed_surface_type, template<typename> class _integration_policy>
struct immersed_structure_traits<ParticleMarkers<_particle_system_type, _immersed_surface_type, _integration_policy> >
{
    typedef typename immersed_structure_traits<_immersed_surface_type>::ode_rhs_type ode_rhs_type;
    typedef typename _particle_system_type::value_type                  value_type;
};

#endif
