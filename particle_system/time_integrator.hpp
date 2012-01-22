#ifndef TIME_INTEGRATOR_HPP
#define TIME_INTEGRATOR_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    TimeIntegrator
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

template<class T >
class surface_traits;
template<typename Derived>
class TimeIntegrator
{
    protected:
        typedef surface_traits<Derived>::time_integrator time_integrator;
        typedef surface_traits<Derived>::value_type value_type;

    public:
      
      TimeIntegrator() {}

        void integrate(value_type time, value_type timestep)
        {
            logger.startTimer("sdcTimeStep");
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

//         inline size_t ode_size() { return derived().data_size(); }     
};



#endif
