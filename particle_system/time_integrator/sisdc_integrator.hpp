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

#include<iterator>
#include "math/ode_solver/sdc/integrator/clenshaw_curtis.hpp"
#include "math/ode_solver/sdc/semi_implicit_sdc.hpp"

template<typename T> struct immersed_structure_traits;

template<typename boundary_type, int sdc_nodes = 5, int sdc_corrections = 4>
class SISDCIntegrator
{
    protected:
        typedef SISDCIntegrator<boundary_type,sdc_nodes,sdc_corrections>               self_type;
        typedef typename immersed_structure_traits<boundary_type>::value_type          value_type;
        typedef SDCSpectralIntegrator<value_type, 0,sdc_nodes,2,sdc_nodes>             spectral_integrator_type;

    private:
        size_t m_ode_size;
        SemiImplicitSDC<value_type,self_type,spectral_integrator_type,sdc_nodes,sdc_corrections> m_sdc;

    public:

        SISDCIntegrator(size_t ode_size) : m_ode_size(ode_size), m_sdc(*this) {}

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
            m_sdc.init(time,positions);
            m_sdc.predictor(time,timestep);
            m_sdc.corrector(time,timestep);
            std::copy(m_sdc.X(0),m_sdc.X(0)+m_ode_size,positions);
            std::transform(m_sdc.Fi(0),m_sdc.Fi(0)+m_ode_size,m_sdc.Fe(0),velocities,std::plus<value_type>());
        }

        void Explicit(value_type time, const value_type *x, value_type *v)
        {
            std::copy(x,x+m_ode_size,boundary().positions());
            boundary().update_forces(time);
            boundary().fluid_solver().Explicit(x,v,boundary().forces());
        }
        
        void Implicit(value_type time, const value_type *x, value_type *v)
        {
            std::copy(x,x+m_ode_size,boundary().positions());
            boundary().update_forces(time);
            boundary().fluid_solver().Implicit(x,v,boundary().forces());
        }
        
        size_t ode_size() { return m_ode_size; }
};

#endif
