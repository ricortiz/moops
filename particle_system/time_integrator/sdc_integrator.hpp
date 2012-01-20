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

#include "math/ode_solver/sdc/integrator/spectral_integrator.hpp"
#include "math/ode_solver/sdc/explicit_sdc.hpp"

template<typename T> struct immersed_structure_traits;

template<typename boundary_type>
class SDCIntegrator
{
    protected:
        typedef SDCIntegrator<boundary_type>                                      self_type;
        typedef typename immersed_structure_traits<boundary_type>::value_type     value_type;
        typedef typename immersed_structure_traits<boundary_type>::ode_rhs_type   ode_rhs_type;
        typedef Integrator<value_type,gauss_lobatto>                              spectral_integrator_type;
        typedef ExplicitSDC<value_type,boundary_type,spectral_integrator_type,9> explicit_sdc_type;

    private:
        explicit_sdc_type  m_sdc;
        ode_rhs_type m_rhs;

    public:

        SDCIntegrator() : m_sdc(derived())
        {
            m_sdc.init(derived().positions(), derived().velocities());
            m_rhs.init(ode_size()/3);
        }

        inline boundary_type &derived()
        {
            return *static_cast<boundary_type*>(this);
        }

        template<typename value_type>
        void integrate(value_type timestep)
        {
            logger.startTimer("sdcTimeStep");
            value_type time = derived().time();
            m_sdc.predictor(time, timestep);
            m_sdc.corrector(time, timestep);
            m_sdc.update();
            logger.stopTimer("sdcTimeStep");
        }

        void operator()(value_type time, const value_type *x, value_type *v)
        {
            logger.startTimer("odeRHSeval");
            derived().updateForces(time);
            m_rhs(x, v, derived().forces());
            logger.stopTimer("odeRHSeval");
        }

        inline size_t ode_size() { return derived().data_size(); }        
        inline ode_rhs_type &ode_rhs() { return m_rhs; }
};



#endif
