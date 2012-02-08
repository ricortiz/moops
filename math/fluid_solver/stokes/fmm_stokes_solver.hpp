#ifndef FMM_STOKES_SOLVER
#define FMM_STOKES_SOLVER
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    FMMStokesSolver
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name FMMStokesSolver - Wrapper for the Fast Multipole Method based Stokes solver
/// @section Description FMMStokesSolver wraps the FMM method so the simulator can use it.
/// @section See also MultipoleTaylor

#include "math/fluid_solver/stokes/fmm/octree/octree.hpp"
#include "math/fluid_solver/stokes/fmm/multipole_taylor.hpp"
#include "particle_system/particle.hpp"

template<typename value_type, size_t max_particles, size_t precision>
class FmmStokesSolver
{
    protected:
        typedef ParticleWrapper<value_type> particle_type;
        typedef Octree<value_type,max_particles,particle_type> tree_type;
        typedef MultipoleTaylor<value_type,precision> fmm_type;
        
    private:
        value_type m_delta;
        size_t m_num_particles;
        
    public:
        FmmStokesSolver(size_t num_particles) : m_num_particles(num_particles) {}
        
        inline void operator()(value_type, value_type *x, value_type *v, value_type *, value_type *f)
        {
            tree_type tree(x,v,f,m_num_particles);
            tree.init();
            fmm_type fmm(tree.boxes().size());
            for (int i = 0; i < 4; ++i)
            {
                fmm.compute_far_field_expansions(tree.root(), i);
                fmm.compute_local_expansions(tree.root(), i);
            }

            fmm.compute_velocity(tree.root(), m_delta);
        }
        
        void setDelta(value_type delta) { m_delta = delta; }
};


#endif

