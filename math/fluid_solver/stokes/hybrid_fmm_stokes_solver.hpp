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
        value_type m_extents[2];
        size_t m_num_particles;
        GPU_Velocities *m_gpu_velocity;
        Particle *m_particles;
        std::vector<Field> m_fields;
        std::vector<Potential> m_potentials;

    public:
        HybridFmmStokesSolver(size_t num_particles) 
	: 
            m_num_particles(num_particles), 
            m_gpu_velocity(new GPU_Velocities[3*num_particles + 2 * MEMORY_ALIGNMENT/sizeof(GPU_Velocities)]),
            m_particles(new Particle[3*num_particles + 2 * MEMORY_ALIGNMENT/sizeof(Particle)]),
            m_fields(num_particles),
            m_potentials(num_particles)
        {
            octree.maxBodiesPerNode = 50;
            octree.numParticles = num_particles;
            int precision = 6;
            octree.fields = &m_fields[0];
	    octree.potentials = &m_potentials[0];
            octree.GPU_Veloc = (GPU_Velocities*) ALIGN_UP(m_gpu_velocity,MEMORY_ALIGNMENT);
            octree.bodies = (Particle *) ALIGN_UP( m_particles, MEMORY_ALIGNMENT );
            float *start = &(octree.GPU_Veloc[0].x);
            std::fill(start,start+3*num_particles,0.0);
            start = &octree.bodies[0].position[0];
            std::fill(start,start+6*num_particles,0.0);
            m_extents[0] = -10; m_extents[1] = 100;
            logger.startTimer("CreateOctree");
            CreateOctree(num_particles, precision, m_extents[1], m_extents[0]);
            logger.stopTimer("CreateOctree");
        }
        ~HybridFmmStokesSolver()
        {
            delete [] m_gpu_velocity;
            delete [] m_particles;            
        }
        inline void operator() ( value_type t, const value_type *x, value_type *v, const  value_type *f)
        {
            
        }
        void setDomain(value_type domain[2][3])
        {
            value_type tmp[3] = {domain[0][0]+domain[1][0],domain[0][1]+domain[1][0],domain[0][2]+domain[1][0]};
            m_extents[0] = 0, m_extents[1] = 0;
            for( int i = 0; i < 3; ++i)
            {
                if(m_extents[0] < tmp[i])
                    m_extents[0] = tmp[i];
                if(m_extents[1] > tmp[i])
                    m_extents[1] = tmp[i];
            }
            
        }
        void setDelta(value_type delta) { m_delta = delta; }
};


#endif

