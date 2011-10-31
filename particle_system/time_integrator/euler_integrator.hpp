#ifndef EULER_INTEGRATOR_HPP
#define EULER_INTEGRATOR_HPP
/*=========================================================================

Program:   Modular Object Oriented Particle Simulator
Module:    EulerIntegrator

Copyright (c) Ricardo Ortiz
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.

=========================================================================*/

/// @name EulerIntegrator - Time integrator based on forward or backward Euler
/// @section Description The EulersIntegrator class update positions of the 
///			 particles using a simple Euler scheme.  Its a superclass of
///			 the elastic boundary class
/// @section See Also ElasticBoundary ParticleMarkers



template<typename T> struct immersed_structure_traits;

template<typename boundary_type>
class EulerIntegrator
{
    protected:
        typedef typename immersed_structure_traits<boundary_type>::value_type          value_type;
        typedef typename immersed_structure_traits<boundary_type>::particle_integrator_type integrator_type;

    private:
        integrator_type euler_integrator;

    public:

        inline boundary_type &boundary()
        {
            return *static_cast<boundary_type*>(this);
        }

        template<typename value_type>
        inline void integrate(value_type timestep)
        {
            value_type time = boundary().time();
            euler_integrator(boundary(),time,timestep);
        }

};



#endif
