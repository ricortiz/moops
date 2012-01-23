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


template<typename Derived>
class TimeIntegrator
{
    protected:
        typedef surface_traits<Derived>::value_type      value_type;

    public:
        inline Derived &derived()
        {
            return *static_cast<Derived*>(this);
        }
        
        inline void integrate(value_type time, value_type timestep)
        {
            logger.startTimer("sdcTimeStep");
	    Derived().eval(time,timestep);
            logger.stopTimer("sdcTimeStep");
        }

        inline void operator()(value_type time, const value_type *x, value_type *v)
        {
            logger.startTimer("odeRHSeval");
	    Derived().eval(time,x,v);
            logger.stopTimer("odeRHSeval");
        }

//         inline size_t ode_size() { return derived().data_size(); }
};



#endif
