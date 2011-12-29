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

template<typename boundary_type>
class SDCIntegrator
{
    protected:
        typedef SDCIntegrator<boundary_type>                       self_type;
        typedef typename immersed_structure_traits<boundary_type>::value_type              value_type;
        typedef typename immersed_structure_traits<boundary_type>::fluid_solver_type       fluid_solver_type;
        typedef SDCSpectralIntegrator<value_type, 0, 5, 2, 5>                  spectral_integrator_type;
        typedef ExplicitSDC<value_type, boundary_type, spectral_integrator_type, 5, 4>  explicit_sdc_type;

    private:
        size_t    m_ode_size;
        explicit_sdc_type  m_sdc;
	fluid_solver_type &m_fluid_solver;

    public:

        SDCIntegrator(size_t ode_size) : m_ode_size(ode_size), m_sdc(derived()), m_fluid_solver(derived().fluid_solver())
        {
            m_sdc.init(derived().positions(), derived().velocities());
        }

        inline boundary_type &derived()
        {
            return *static_cast<boundary_type*>(this);
        }

        template<typename value_type>
        void integrate(value_type timestep)
        {
            value_type time = derived().time();
            m_sdc.predictor(time, timestep);
            m_sdc.corrector(time, timestep);
        }

        void operator()(value_type time, value_type *x, value_type *v)
        {
            std::copy(x, x + m_ode_size, derived().positions());
            derived().update_forces(time);
            m_fluid_solver(x, v, derived().forces());
        }

        size_t ode_size() { return m_ode_size; }

};



#endif
