#ifndef HYBRID_FMM_STOKES_SOLVER_HPP
#define HYBRID_FMM_STOKES_SOLVER_HPP

#include <vector>

extern "C"
{
#include "math/fluid_solver/stokes/fmm/hybrid/hybridfmm.h"
}

template<typename value_type>
class HybridFmmStokesSolver
{
    private:
        size_t m_num_particles;
        float *m_gpu_velocity;
        Particle *m_particles;
        std::vector<Field> m_fields;
        std::vector<Potential> m_potentials;
        value_type m_delta;
        value_type m_extents[2];

    public:
        HybridFmmStokesSolver(size_t num_particles)
                :
                m_num_particles(num_particles),
                m_gpu_velocity(new float[3*num_particles + 2 * MEMORY_ALIGNMENT/sizeof(float)]),
                m_particles(new Particle[3*num_particles + 2 * MEMORY_ALIGNMENT/sizeof(Particle)]),
                m_fields(num_particles),
                m_potentials(num_particles)
        {
            octree.maxBodiesPerNode = 50;
            octree.numParticles = num_particles;
            octree.fields = &m_fields[0];
            octree.potentials = &m_potentials[0];
            octree.GPU_Veloc = (float*) ALIGN_UP(m_gpu_velocity, MEMORY_ALIGNMENT);
            octree.bodies = (Particle *) ALIGN_UP( m_particles, MEMORY_ALIGNMENT );
            float *start = octree.GPU_Veloc;
            std::fill(start, start + 3*num_particles, 0.0);
            start = &octree.bodies[0].position[0];
            std::fill(start, start + 6*num_particles, 0.0);

        }
        ~HybridFmmStokesSolver()
        {
            delete [] m_gpu_velocity;
            delete [] m_particles;
        }
        inline void operator() ( value_type, const value_type *x, value_type *v, const  value_type *f)
        {
            initData(x,v,f);
            logger.startTimer("CreateOctree");
            CreateOctree(m_num_particles, 6, 256, 0);
            logger.stopTimer("CreateOctree");
//             #pragma omp parallel
            {
                #pragma omp single
                {
                    logger.startTimer("gpuVelocitiesEval");
                    gpuVelocitiesEval(octree.numParticles, octree.numLeafDInodes,
                    octree.total_interaction_pairs, (float *) octree.bodies,
                    (float *) octree.GPU_Veloc, octree.target_list,
                                      octree.number_IP, octree.interaction_pairs,m_delta);
                    logger.stopTimer("gpuVelocitiesEval");
                }
                #pragma omp single
                {                    
                    UpSweep(octree.root);
                    DownSweep(octree.root);
                }
            }
            #pragma omp barrier
            copyVelocities(v);
        }

        void initData(const value_type *x, value_type *v, const value_type *f)
        {
            for (size_t p = 0, idx = 0; p < m_num_particles; ++p, idx += 3)
            {
                octree.bodies[p].position[0] = x[idx];
                octree.bodies[p].position[1] = x[idx+1];
                octree.bodies[p].position[2] = x[idx+2];
                octree.bodies[p].force[0] = f[idx];
                octree.bodies[p].force[1] = f[idx+1];
                octree.bodies[p].force[2] = f[idx+2];
            }
            octree.CPU_Veloc = v;
        }

        void copyVelocities(value_type *v)
        {
            for (size_t p = 0, idx = 0; p < m_num_particles; ++p, idx += 3)
            {
                v[idx] += octree.GPU_Veloc[idx];
                v[idx+1] += octree.GPU_Veloc[idx+1];
                v[idx+2] += octree.GPU_Veloc[idx+2];
            }
        }
        void setDomain(value_type domain[2][3])
        {
            value_type tmp[3] = {domain[0][0] + domain[1][0], domain[0][1] + domain[1][0], domain[0][2] + domain[1][0]};
            m_extents[0] = 0, m_extents[1] = 0;
            for ( int i = 0; i < 3; ++i)
            {
                if (m_extents[0] < tmp[i])
                    m_extents[0] = tmp[i];
                if (m_extents[1] > tmp[i])
                    m_extents[1] = tmp[i];
            }

        }
        void setDelta(value_type delta) { m_delta = delta; }
};


#endif

