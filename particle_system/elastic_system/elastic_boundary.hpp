#ifndef ELASTIC_BOUNDARY_HPP
#define ELASTIC_BOUNDARY_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    ElasticBoundary
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name ElasticBoundary
/// @section Description
/// @section See also

#include <vector>
#include "geometry/surface.hpp"

template < typename spring_system_type,
typename ode_rhs_type,
template<typename> class integration_policy >
class ElasticBoundary :
            public spring_system_type,
            public integration_policy<ElasticBoundary<spring_system_type, ode_rhs_type, integration_policy> >
{
    public:
        typedef ElasticBoundary<spring_system_type, ode_rhs_type, integration_policy> self_type;
        typedef typename spring_system_type::value_type    value_type;
        typedef typename spring_system_type::spring_type   spring_type;
        typedef typename spring_system_type::particle_type particle_type;
        typedef integration_policy<self_type>              time_integrator_type;
	
    private:
        size_t    m_ode_size;

    public:
        ElasticBoundary(size_t ode_size) : m_ode_size(ode_size), spring_system_type(ode_size) {  }
        ~ElasticBoundary() {}

        void run(value_type timestep)
        {
            this->integrate(timestep);
            this->time() += timestep;
        }
};

template< typename _spring_system_type, typename _ode_rhs_type, template<typename> class _integration_policy>
struct immersed_structure_traits<ElasticBoundary<_spring_system_type, _ode_rhs_type, _integration_policy> >
{
    typedef _ode_rhs_type ode_rhs_type;
    typedef typename _spring_system_type::value_type    value_type;
};

#endif
