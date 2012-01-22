#ifndef EXPLICIT_SDC_TIME_INTEGRATOR_HPP
#define EXPLICIT_SDC_TIME_INTEGRATOR_HPP
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

#include "sdc/integrator/spectral_integrator.hpp"
#include "sdc/explicit_sdc.hpp"

template<typename value_type, typename rhs_type>
class SdcTimeIntegrator
{
    protected:
        typedef Integrator<value_type,gauss_lobatto>                              spectral_integrator_type;
        typedef ExplicitSDC<value_type,SdcTimeIntegrator,spectral_integrator_type,9> explicit_sdc_type;

    private:
        rhs_type  m_rhs;
        explicit_sdc_type  m_sdc;

    public:
      
      SdcTimeIntegrator(rhs_type &rhs) : m_rhs(rhs),m_sdc(*this)  {}

	void init(value_type *positions, value_type *velocities)
	{
	  m_sdc.init(positions,velocities);
	}
	
        void integrate(value_type time, value_type timestep)
        {
            logger.startTimer("sdcTimeStep");
            m_sdc.predictor(time, timestep);
            m_sdc.corrector(time, timestep);
            m_sdc.update();
            logger.stopTimer("sdcTimeStep");
        }

//         template<typename value_type>
//         void operator()(value_type time, const value_type *x, value_type *v)
//         {
//             logger.startTimer("odeRHSeval");
//             derived().updateForces(time);
//             m_rhs(x, v, derived().forces());
//             logger.stopTimer("odeRHSeval");
//         }

//         inline size_t ode_size() { return derived().data_size(); }     
};



#endif
