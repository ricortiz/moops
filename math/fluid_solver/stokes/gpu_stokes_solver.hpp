#ifndef GPU_STOKES_SOLVER
#define GPU_STOKES_SOLVER
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    CudaStokesSolver
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name CudaStokesSolver - Wrapper for the direct cuda-based stokes solver.
/// @section Description Wraps the cuda kernels that solve Stokes equations so 
///                      it can be used by the simulator.
/// @section See also cuda

#include "math/fluid_solver/stokes/nbody_gpu/gpu_compute_velocity.hpp"

template<typename particle_system_type>
class GpuStokesSolver
{
    protected:
        typedef typename particle_system_type::value_type value_type;
        
    private:

        value_type m_delta;
        value_type m_num_particles;

    public:

        GpuStokesSolver(value_type num_particles) : m_num_particles(num_particles) {}
        
        inline void operator()(value_type *x, value_type *v, value_type *f)
        {
            operator()(x,v,x,f,m_num_particles,m_num_particles);
        }

        inline void operator()(value_type *x, value_type *v, value_type *y, value_type *f, size_t num_sources, size_t num_targets)
        {
            ComputeStokeslets(x,v,y,f,m_delta,num_sources,num_targets);
        }

        value_type const &delta() const { return m_delta; }
        value_type &delta() { return m_delta; }
        
};


#endif

