#ifndef PARTICLE_MARKERS_HPP
#define PARTICLE_MARKERS_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    ParticleMarkers
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name ParticleMarkers - Stores data for freely moving points in the fluid.
/// @section Description The ParticleMarkers class is used in the case you want
///     to move points in the fluid that have no constraints
///     (ie, there is no springs connecting any of the points)


template < typename particle_system_type,
typename immersed_surface_type,
template<typename> class integration_policy >
class ParticleMarkers : public particle_system_type, public integration_policy<ParticleMarkers<particle_system_type, immersed_surface_type, integration_policy> >
{
    public:
        typedef typename immersed_structure_traits<immersed_surface_type>::fluid_solver_type fluid_solver_type;
        typedef typename particle_system_type::value_type    value_type;
        typedef typename particle_system_type::particle_type particle_type;
        typedef integration_policy<ParticleMarkers<particle_system_type,immersed_surface_type,integration_policy> > time_integrator_type;

    private:
        immersed_surface_type            &m_surface;
        fluid_solver_type                &m_fluid_solver;

    public:

        ParticleMarkers(immersed_surface_type &surface, size_t ode_size) : particle_system_type(ode_size), time_integrator_type(ode_size), m_surface(surface), m_fluid_solver(surface.fluid_solver()) {  }
        ~ParticleMarkers() {}

        inline fluid_solver_type &fluid_solver() { return m_fluid_solver; }
        inline fluid_solver_type const &fluid_solver() const { return m_fluid_solver; }

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
            m_fluid_solver(x, v, y, f, num_sources, num_targets);
        }
};

template< typename _particle_system_type, typename _immersed_surface_type, template<typename> class _integration_policy>
struct immersed_structure_traits<ParticleMarkers<_particle_system_type, _immersed_surface_type, _integration_policy> >
{
    typedef typename immersed_structure_traits<_immersed_surface_type>::fluid_solver_type fluid_solver_type;
    typedef typename _particle_system_type::value_type                  value_type;
};

#endif
