#ifndef PARTICLE_MARKERS_HPP
#define PARTICLE_MARKERS_HPP

/// Author: Ricardo Ortiz <ricardo.ortiz@tulane.edu>, (C) 2008
/// $Id: stokes_flow.h 19 2008-12-17 23:41:26Z rortiz $

/**
 * Immersed volume that interact with immersed boundary
 *
 * Takes a geometry and builds a particle system from it.
 */
template<typename particle_system_type,
         typename immersed_surface_type,
         typename fluid_solver_type,
         template<typename> class integration_policy>
         class ParticleMarkers : public particle_system_type, public integration_policy<ParticleMarkers<particle_system_type,immersed_surface_type,fluid_solver_type,integration_policy> >
{ 
    public:
        typedef typename particle_system_type::value_type    value_type;
        typedef typename particle_system_type::particle_type particle_type;

    private:
        fluid_solver_type                &m_fluid_solver;
        immersed_surface_type            *m_surface;

    public:
        
        ParticleMarkers(immersed_surface_type *surface) : m_surface(surface), m_fluid_solver(surface->fluid_solver()) {  }
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
            value_type *f = m_surface->forces();
            value_type *y = m_surface->positions();
            size_t num_targets = this->particles_size();
            size_t num_sources = m_surface->particles_size();
            m_fluid_solver(x,v,y,f,num_sources,num_targets);
        }
};

template< typename _particle_system_type, typename _immersed_surface_type, typename _fluid_solver_type, template<typename> class _integration_policy>
struct immersed_structure_traits<ParticleMarkers<_particle_system_type,_immersed_surface_type,_fluid_solver_type,_integration_policy> >
{
    typedef typename _particle_system_type::value_type                  value_type;
    typedef typename _particle_system_type::particle_integrator_type    particle_integrator_type;
};

#endif
