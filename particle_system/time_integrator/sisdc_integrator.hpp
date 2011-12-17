#ifndef SEMI_IMPLICIT_SDC_INTEGRATOR_HPP
#define SEMI_IMPLICIT_SDC_INTEGRATOR_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    SISDCIntegrator
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name SISDCIntegrator -- Wrapper for the Spectral Deffered Correction time integrator
/// @section Description SISDCIntegrator wraps the integrator so it can be used by the simulator
///                      It contains the main method that drives the simulation: void integrate()
/// @section See also ExplicitSDC SemiImplicitSDC

#include "math/ode_solver/sdc/sdc_base.hpp"

template<typename T> struct immersed_structure_traits;

template<typename boundary_type>
class SISDCIntegrator
{
    protected:
        typedef typename immersed_structure_traits<boundary_type>::value_type          value_type;
        typedef typename immersed_structure_traits<boundary_type>::particle_integrator_type sdc_integrator_type;

    private:
        enum
        {
            m_ode_size = sdc_traits<sdc_integrator_type>::ode_size
        };
        sdc_integrator_type m_sdc;

    public:

        SISDCIntegrator() {}

        inline boundary_type &boundary()
        {
            return *static_cast<boundary_type*>(this);
        }

        template<typename value_type>
        void integrate(value_type timestep)
        {
            value_type *positions = boundary().positions();
            value_type *velocities = boundary().velocities();
            value_type time = boundary().time();
            std::vector<value_type> Vi(m_ode_size,0.0);
            std::vector<value_type> Ve(m_ode_size,0.0);
            boundary().Implicit(time,positions,&Vi[0]);
            boundary().Explicit(time,positions,&Ve[0]);
            m_sdc.init(positions,&Vi[0],&Ve[0]);
            m_sdc.predictor(boundary(),time,timestep);
            m_sdc.corrector(boundary(),time,timestep);
            std::copy(m_sdc.X(0),m_sdc.X(0)+m_ode_size,positions);
            std::transform(m_sdc.Fi(0),m_sdc.Fi(0)+m_ode_size,m_sdc.Fe(0),velocities,std::plus<value_type>());
        }
};



#endif
