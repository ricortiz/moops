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
#include "math/ode_solver/sdc/integrator/spectral_integrator.hpp"
#include "math/ode_solver/sdc/semi_implicit_sdc.hpp"

template<typename T> struct immersed_structure_traits;

template<typename boundary_type>
class SISDCIntegrator
{
    protected:
        typedef SISDCIntegrator<boundary_type>               self_type;
        typedef typename immersed_structure_traits<boundary_type>::value_type              value_type;
        typedef typename immersed_structure_traits<boundary_type>::ode_rhs_type       ode_rhs_type;
        typedef Integrator<value_type,gauss_lobatto>                  spectral_integrator_type;
	typedef SemiImplicitSDC<value_type, boundary_type, spectral_integrator_type, 9> semi_implicit_sdc_type;
	
    private:
        size_t m_ode_size;
        semi_implicit_sdc_type m_sisdc;
        ode_rhs_type m_rhs;
        std::vector<value_type> m_Ve0;
        std::vector<value_type> m_Vi0;

    public:
        SISDCIntegrator(size_t ode_size) : m_ode_size(ode_size), m_sisdc(derived()), m_rhs(ode_size/3), m_Vi0(ode_size, 0.0), m_Ve0(ode_size, 0.0)
        {
            m_sisdc.init(derived().positions(), &m_Vi0[0], &m_Ve0[0]);
            m_rhs.initMaps(derived());
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
            m_sisdc.predictor(time, timestep);
            m_sisdc.corrector(time, timestep);
            m_sisdc.update();
            std::transform(m_Vi0.begin(), m_Vi0.end(), m_Ve0.begin(), derived().velocities(), std::plus<value_type>());
            logger.stopTimer("sdcTimeStep");
        }

        void Explicit(value_type time, const value_type *x, value_type *v)
        {
            logger.startTimer("odeRHSExplicit");
            derived().updateForces(time,x);
            m_rhs.Explicit(x, v, derived().forces());
            logger.stopTimer("odeRHSExplicit");
        }

        void Implicit(value_type time, const value_type *x, value_type *v)
        {
            logger.startTimer("odeRHSImplicit");
            derived().updateForces(time,x);
            m_rhs.Implicit(x, v, derived().forces());
            logger.stopTimer("odeRHSImplicit");
        }

        size_t ode_size() { return m_ode_size; }

        ode_rhs_type &ode_rhs() { return m_rhs; }
};

#endif
