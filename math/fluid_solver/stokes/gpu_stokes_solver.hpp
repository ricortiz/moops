#ifndef GPU_STOKES_SOLVER
#define GPU_STOKES_SOLVER

#include "math/fluid_solver/stokes/nbody_gpu/gpu_compute_velocity.hpp"

template<typename value_type>
class GpuStokesSolver
{        
    private:
        value_type m_delta;
        size_t m_num_particles;

    public:
        GpuStokesSolver(size_t num_particles) : m_num_particles(num_particles) {}
        
        inline void operator()(value_type *x, value_type *v, value_type *f)
        {
            operator()(x,v,x,f,m_num_particles,m_num_particles);
        }

        inline void operator()(value_type *x, value_type *v, value_type *y, value_type *f, size_t num_sources, size_t num_targets)
        {
            std::fill(v,v+num_targets*3,0.0);
            ComputeStokeslets(x,v,y,f,m_delta,num_sources,num_targets);
        }

        value_type const &setDelta(value_type delta) const { return m_delta = delta; }
        
};


#endif

