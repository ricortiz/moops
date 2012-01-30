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
        typedef typename Traits<Derived>::time_integrator_type time_integrator_type;

    private:
        time_integrator_type time_integrator;

    public:
        TimeIntegrator(size_t ode_size) : time_integrator(ode_size) {}
        
        inline Derived &derived()
        {
            return *static_cast<Derived*>(this);
        }

        template<typename value_type>
        inline void integrate(value_type t, value_type timestep)
        {
            time_integrator(derived(),t,derived().positions(),derived().velocities(),timestep);
        }
        
};



#endif
