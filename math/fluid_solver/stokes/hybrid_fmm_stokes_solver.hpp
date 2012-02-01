#ifndef HYBRID_FMM_STOKES_SOLVER_HPP
#define HYBRID_FMM_STOKES_SOLVER_HPP

#include <vector>

extern "C" {
#include "math/fluid_solver/stokes/fmm/hybrid/hybridfmm.h"
}

template<typename value_type>
class HybridFmmStokesSolver
{
    private:
        value_type m_delta;
        map_type m_map_explicit;
        map_type m_map_implicit;
        size_t m_num_particles;
        std::vector<float> m_gpu_velocity;
        std::vector<Potential> m_fields;
        std::vector<Field> m_potentials;

    public:
        HybridFmmStokesSolver(size_t num_particles) 
	: 
	m_num_particles(num_particles), 
	m_gpu_velocity(3*num_particles + 1 << 11),
	m_fields(num_particles),
	m_potentials(num_particles)	
        {
            octree.maxBodiesPerNode = 50;
            octree.numParticles = num_particles;
            int precision = 6;
            octree.fields = &m_fields[0];
	    octree.potentials = &m_potentials[0];
	    octree.GPU_Veloc = &m_gpu_velocity[0];
            
        }
        
        inline void operator() ( value_type t, const value_type *x, value_type *v, const  value_type *f)
        {
	  octree.CPU_Veloc = v;
        }



        void setDelta(value_type delta) { m_delta = delta; }
};


#endif

