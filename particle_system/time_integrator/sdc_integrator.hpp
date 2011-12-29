#ifndef EXPLICIT_SDC_INTEGRATOR_HPP
#define EXPLICIT_SDC_INTEGRATOR_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    SDCIntegrator
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name SDCIntegrator -- Wrapper for the Spectral Deffered Correction time integrator
/// @section Description SDCIntegrator wraps the integrator so it can be used by the simulator
///                      It contains the main method that drives the simulation: void integrate()
/// @section See also ExplicitSDC SemiImplicitSDC

#include "math/ode_solver/sdc/integrator/clenshaw_curtis.hpp"
#include "math/ode_solver/sdc/explicit_sdc.hpp"

template<typename T> struct immersed_structure_traits;

template<typename boundary_type, int sdc_nodes = 5, int sdc_corrections = 4 >
class SDCIntegrator
{
    protected:
        typedef SDCIntegrator<boundary_type,sdc_nodes,sdc_corrections>                 self_type;
        typedef typename immersed_structure_traits<boundary_type>::value_type          value_type;        
        typedef SDCSpectralIntegrator<value_type, 0,sdc_nodes,2,sdc_nodes>             spectral_integrator_type;

    private:
        size_t m_ode_size;
        ExplicitSDC<value_type,self_type,spectral_integrator_type,sdc_nodes,sdc_corrections> m_sdc;

    public:

        SDCIntegrator(size_t ode_size) : m_ode_size(ode_size), m_sdc(*this) {}
        
        inline boundary_type &derived()
        {
            return *static_cast<boundary_type*>(this);
        }

        template<typename value_type>
        void integrate(value_type timestep)
        {
            value_type *positions = derived().positions();
            value_type *velocities = derived().velocities();
            value_type time = derived().time();
            m_sdc.init(time,positions);
            m_sdc.predictor(time,timestep);
            m_sdc.corrector(time,timestep);
            std::copy(m_sdc.X(0),m_sdc.X(0)+m_ode_size,positions);
            std::copy(m_sdc.F(0),m_sdc.F(0)+m_ode_size,velocities);
        }

        void operator()(value_type time, value_type *x, value_type *v)
        {
            std::copy(x,x+m_ode_size,derived().positions());
            derived().update_forces(time);
            derived().fluid_solver()(x,v,derived().forces());
        }

        size_t ode_size() { return m_ode_size; }

};



#endif
