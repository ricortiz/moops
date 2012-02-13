#ifndef HYBRID_FMM_STOKES_SOLVER_HPP
#define HYBRID_FMM_STOKES_SOLVER_HPP

#include <vector>

extern "C"
{
#include "math/fluid_solver/stokes/fmm/hybrid/hybridfmm.h"
}

template < typename value_type, int precision = 9 >
class HybridFmmStokesSolver
{
    private:
        size_t m_num_particles;
        std::vector<float> m_gpu_velocity;
        std::vector<Particle> m_particles;
        std::vector<Field> m_fields;
        std::vector<Potential> m_potentials;
        value_type m_delta;
        value_type m_extents[2];
        bool m_initialized;

    public:
        HybridFmmStokesSolver(size_t num_particles)
                :
                m_num_particles(num_particles),
                m_gpu_velocity(3*num_particles),
                m_particles(num_particles),
                m_fields(num_particles),
                m_potentials(num_particles),
                m_initialized(false)

        {
            octree.maxBodiesPerNode = 2;
            octree.numParticles = num_particles;
            octree.fields = &m_fields[0];
            octree.potentials = &m_potentials[0];
            octree.GPU_Veloc = &m_gpu_velocity[0];
            octree.bodies = &m_particles[0];
        }
        inline void operator()(value_type, value_type *x, value_type *v, value_type *f)
        {
            initOctree(x, v, f);
// #pragma omp parallel
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
                    std::cout << "fmm_velocities = [";std::copy(v, v + 3*m_num_particles, std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "]" << std::endl;
                }
#pragma omp barrier
#pragma omp single
                {
                    gpuGetVelocities();
                    std::cout << "gpu_velocities = [";std::copy(m_gpu_velocity.begin(), m_gpu_velocity.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "]" << std::endl;
                }
#pragma omp barrier
// #pragma omp single
//                 {
//                     //ReSort(octree.root, octree.rootInfo);
//                     /*if(ReSort(octree.root, octree.rootInfo))
//                      *       RebuildTree(number_particles, precision, maximum_extent, minimum_extent);
//                      * */
//                 }
            }
            logger.startTimer("copyVelocities");
            copyVelocities(v);
            logger.stopTimer("copyVelocities");
            std::cout << "hv_velocities = [";std::copy(v, v + 3*m_num_particles, std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "]" << std::endl;
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
                octree.bodies[p].index = p;
            }
        }

        inline void initOctree(const value_type *x, value_type *v, const value_type *f)
        {
            logger.startTimer("initData");
            initData(x, f);
            print_particle_index();
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
                updateOctree();
                logger.stopTimer("UpdateOctree");
            }
            print_particle_index();
        }

        void copyVelocities(value_type *v)
        {
            std::transform(m_gpu_velocity.begin(), m_gpu_velocity.end(), v, m_gpu_velocity.end(), std::plus<value_type>());
            size_t p = 0, idx = 0;
// #pragma omp parallel for shared(v) private(p) firstprivate(idx)
            for (p = 0; p < m_num_particles; ++p)
            {
                size_t idx_p = 3 * octree.bodies[p].index;
                v[idx_p] = m_gpu_velocity[idx];
                v[idx_p+1] = m_gpu_velocity[idx+1];
                v[idx_p+2] = m_gpu_velocity[idx+2];
                idx += 3;
            }
            std::fill(m_gpu_velocity.begin(), m_gpu_velocity.end(), 0.0);
        }

        
        void setDelta(value_type delta) { m_delta = delta; }

        void print_particle_index()
        {
            std::cout << "particle = [";
            for (size_t i = 0; i < m_num_particles; ++i)
                std::cout << octree.bodies[i].index << " ";
            std::cout << "]" << std::endl;
        }

        
        void allPairs()
        {
            std::vector<value_type> velocity(3*m_num_particles, 0.0);
            for (size_t i = 0; i < m_num_particles; ++i)
            {
                unsigned i_idx = octree.bodies[i].index;
                for (size_t j = 0; j < m_num_particles; ++j)
                    ComputeVelocityP2P(octree.bodies[i_idx].position, &velocity[i_idx], octree.bodies[j].position, octree.bodies[j].force);
            }
            std::cout << "ap_velocities = [";std::copy(velocity.begin(), velocity.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "]" << std::endl;
        }
        
        void ComputeVelocityP2P(value_type *target, value_type *velocity, value_type *source, value_type *force)
        {

            value_type dx[3] = {target[0] - source[0], target[1] - source[1], target[2] - source[2]};

            value_type r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            value_type d2 = m_delta * m_delta;
            value_type R1 = r2 + d2;
            value_type R2 = R1 + d2;
            value_type invR = 1.0 / R1;
            value_type H = sqrt(invR) * invR * 0.039788735772974;

            value_type fdx = force[0] * dx[0] + force[1] * dx[1] + force[2] * dx[2];
            velocity[0] += H * (force[0] * R2 + fdx * dx[0]);
            velocity[1] += H * (force[1] * R2 + fdx * dx[1]);
            velocity[2] += H * (force[2] * R2 + fdx * dx[2]);
        }     
};


#endif


