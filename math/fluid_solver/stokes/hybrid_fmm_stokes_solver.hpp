#ifndef HYBRID_FMM_STOKES_SOLVER_HPP
#define HYBRID_FMM_STOKES_SOLVER_HPP

#include <vector>
#include "cpu_stokes_solver.hpp"

extern "C"
{
#include "math/fluid_solver/stokes/fmm/hybrid/hybridfmm.h"
}

template < typename value_type, int precision = 6 >
class HybridFmmStokesSolver
{
    private:
        size_t m_num_particles;
        std::vector<float> m_gpu_velocity;
        std::vector<Particle> m_particles;
        std::vector<Field> m_fields;
        std::vector<Potential> m_potentials;
        std::vector<unsigned int> m_index;
        value_type m_delta;
        bool m_initialized;

    public:
        HybridFmmStokesSolver(size_t num_particles)
                :
                m_num_particles(num_particles),
                m_gpu_velocity(3*num_particles),
                m_particles(num_particles),
                m_fields(num_particles),
                m_potentials(num_particles),
                m_index(num_particles),
                m_initialized(false)
        {
            octree.maxBodiesPerNode = 2;
            octree.numParticles = num_particles;
            octree.fields = &m_fields[0];
            octree.potentials = &m_potentials[0];
            octree.GPU_Veloc = &m_gpu_velocity[0];
            octree.bodies = &m_particles[0];
	    octree.particle_idx = &m_index[0];
        }
        inline void operator()(value_type, value_type *x, value_type *v, value_type *f)
        {
            initOctree(x,v,f);
	    #pragma omp parallel
            {
		#pragma omp single
                {
                    logger.startTimer("gpuVelocitiesEval");
                    gpuVelocitiesEval(octree.numParticles,
                                      octree.numLeafDInodes,
                                      octree.total_interaction_pairs,
                                      (float*)octree.bodies,
                                      octree.GPU_Veloc,
                                      octree.target_list,
                                      octree.number_IP,
                                      octree.interaction_pairs,
                                      m_delta);
                    logger.stopTimer("gpuVelocitiesEval");
                }
		#pragma omp single
                {
                    logger.startTimer("UpSweep");
                    UpSweep(octree.root);
                    logger.stopTimer("UpSweep");
                    logger.startTimer("DownSweep");
                    DownSweep(octree.root);
                    logger.stopTimer("DownSweep");
//                     std::cout << "fmm_velocities = ";std::copy(v, v + 3*m_num_particles, std::ostream_iterator<value_type>(std::cout, " ")); std::cout << std::endl;
                }
		#pragma omp barrier
		#pragma omp single
                {
                    gpuGetVelocities();
                }
		#pragma omp barrier
            }
//             std::cout << "gpu_velocities = ";std::copy(m_gpu_velocity.begin(), m_gpu_velocity.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << std::endl;
            logger.startTimer("copyVelocities");
            copyVelocities(v);
            logger.stopTimer("copyVelocities");
            std::cout << "hv_velocities = ";std::copy(v, v + 3*m_num_particles, std::ostream_iterator<value_type>(std::cout, " ")); std::cout << std::endl;
        }

        void initData(const value_type *x, const value_type *f)
        {
            for (size_t p = 0, idx = 0; p < m_num_particles; ++p, idx += 3)
            {
                octree.bodies[p].position[0] = x[idx];
                octree.bodies[p].position[1] = x[idx+1];
                octree.bodies[p].position[2] = x[idx+2];
                octree.bodies[p].force[0] = f[idx];
                octree.bodies[p].force[1] = f[idx+1];
                octree.bodies[p].force[2] = f[idx+2];
		octree.particle_idx[p] = p;
            }
        }

        inline void initOctree(const value_type *x, value_type *v, const value_type *f)
        {
            logger.startTimer("initData");
            initData(x, f);
// 	    print_particle_index();
            octree.CPU_Veloc = v;
            logger.stopTimer("initData");
            if (!m_initialized)
            {
                logger.startTimer("CreateOctree");
                CreateOctree(m_num_particles, precision, 255.9999, 0);
                logger.stopTimer("CreateOctree");
                m_initialized = true;
            }
            else
            {
                logger.startTimer("UpdateOctree");
//                 if (ReSort(octree.root, octree.rootInfo))
                    RebuildTree(m_num_particles, precision, 255.9999, 0);
                logger.stopTimer("UpdateOctree");
            }
//             print_particle_index();
        }
        void print_particle_index()
        {
            std::cout << "particle_idx = [";std::copy(m_index.begin(), m_index.end(), std::ostream_iterator<unsigned int>(std::cout, " ")); std::cout << "]" << std::endl;
        }
        void copyVelocities(value_type *v)
        {
            std::transform(m_gpu_velocity.begin(), m_gpu_velocity.end(), v, m_gpu_velocity.begin(), std::plus<value_type>());
            size_t p = 0, idx = 0;
            for (p = 0; p < m_num_particles; ++p)
            {
                size_t idx_p = 3 * octree.particle_idx[p];
                v[idx_p] = m_gpu_velocity[idx];
                v[idx_p+1] = m_gpu_velocity[idx+1];
                v[idx_p+2] = m_gpu_velocity[idx+2];
                idx += 3;
            }
	    std::fill(m_gpu_velocity.begin(), m_gpu_velocity.end(), 0.0);
        }

        void allPairs()
        {
            std::vector<value_type> velocity(3*m_num_particles, 0.0);
            for (size_t i = 0; i < m_num_particles; ++i)
            {                
                unsigned i_idx = 3*octree.particle_idx[i];
                for (size_t j = 0; j < m_num_particles; ++j)
                    compute_velocity(octree.bodies[i].position, &velocity[i_idx], octree.bodies[j].position, octree.bodies[j].force,m_delta);
            }
//             std::cout << "ap_velocities = [";std::copy(velocity.begin(), velocity.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "]" << std::endl;
        }
        void setDelta(value_type delta) { m_delta = delta; }
};


#endif


