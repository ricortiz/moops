#ifndef IMMERSED_BOUNDARY_HPP
#define IMMERSED_BOUNDARY_HPP

/// Author: Ricardo Ortiz <ricardo.ortiz@tulane.edu>, (C) 2008
/// $Id: stokes_flow.h 19 2008-12-17 23:41:26Z rortiz $

#include <vector>
#include "geometry/surface.hpp"
/**
 * Immersed boundary that interact with fluid
 *
 * Takes a geometry and builds a particle system from it.
 */
template<typename spring_system_type, typename fluid_solver_type, template<typename> class integration_policy>
class ImmersedBoundary : public spring_system_type, public integration_policy<ImmersedBoundary<spring_system_type,fluid_solver_type,integration_policy> >
{ 
    public:
        typedef typename spring_system_type::value_type    value_type;
        typedef typename spring_system_type::spring_type   spring_type;
        typedef typename spring_system_type::particle_type particle_type;

    private:
        fluid_solver_type                m_fluid_solver;

    public:
        ImmersedBoundary() {  }
        ~ImmersedBoundary() {}
        inline fluid_solver_type &fluid_solver() { return m_fluid_solver; }
        inline fluid_solver_type const &fluid_solver() const { return m_fluid_solver; }
        void run(value_type timestep)
        {
            this->integrate(timestep);
            this->time() += timestep;
        }

        void operator()(value_type time, value_type *x, value_type *v)
        {
            this->update_forces(time);
            value_type *f = this->forces();
            m_fluid_solver(x,v,f);
        }
};

template< typename _spring_system_type, typename _fluid_solver_type, template<typename> class _integration_policy>
struct immersed_structure_traits<ImmersedBoundary<_spring_system_type,_fluid_solver_type,_integration_policy> >
{
    typedef typename _spring_system_type::value_type    value_type;
    typedef typename _spring_system_type::particle_integrator_type    particle_integrator_type;
};

#endif
