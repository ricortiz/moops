#ifndef DIRECT_STOKES_SOLVER_HPP
#define DIRECT_STOKES_SOLVER_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    DirectStokesSolver
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name DirectStokesSolver - Solves the Stokes problem using a direct solver.
/// @section Description DirectStokesSolver solves an n-body problem and updates the
///                      velocities of the particle system.
/// @section See also 
#include "math/fluid_solver/stokes/nbody_cpu/cpu_compute_velocity.hpp"

template<typename particle_system_type>
class DirectStokesSolver
{
    protected:
        typedef typename particle_system_type::value_type value_type;
        
    private:

        value_type m_delta;

        enum
        {
            num_particles = particle_system_type::num_particles
        };


    public:
        inline void operator()(value_type *x, value_type *v, value_type *f)
        {
            operator()(x,v,x,f,num_particles,num_particles);
        }

        inline void operator()(value_type *x, value_type *v, value_type *y, value_type *f, size_t num_sources, size_t num_targets)
        {
            size_t size_targets = 3*num_targets;
            size_t size_sources = 3*num_sources;
            for (size_t i = 0; i < size_targets; i+=3)
                for (size_t j = 0; j < size_sources; j+=3)
                    compute_velocity(&x[i],&v[i],&y[j],&f[j],m_delta);
        }
   
        value_type const &delta() const { return m_delta; }
        value_type &delta() { return m_delta; }
};


#endif

