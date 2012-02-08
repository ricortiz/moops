#ifndef HYBRID_FMM_STOKES_SOLVER_HPP
#define HYBRID_FMM_STOKES_SOLVER_HPP

#include <vector>

extern "C"
{
#include "math/fluid_solver/stokes/fmm/hybrid/hybridfmm.h"
}

template<typename value_type, int precision = 6>
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

    public:
        HybridFmmStokesSolver(size_t num_particles)
                :
                m_num_particles(num_particles),
                m_gpu_velocity(3*num_particles),
                m_particles(num_particles),
                m_fields(num_particles),
                m_potentials(num_particles)
        {
            octree.maxBodiesPerNode = 25;
            octree.numParticles = num_particles;
            octree.fields = &m_fields[0];
            octree.potentials = &m_potentials[0];
            octree.GPU_Veloc = &m_gpu_velocity[0];
            octree.bodies = &m_particles[0];
        }
        inline void operator()(value_type, value_type *x, value_type *v, value_type *f)
        {
            logger.startTimer("initData");
            initData(x, v, f);
            logger.stopTimer("initData");
            logger.startTimer("CreateOctree");
            CreateOctree(m_num_particles, precision, 255.9999, 0);
            logger.stopTimer("CreateOctree");
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
                    DownSweep(octree.root);
                    logger.stopTimer("UpSweep");
//                     std::cout << "fmm_velocities = ";std::copy(v, v + 3*m_num_particles, std::ostream_iterator<value_type>(std::cout, " ")); std::cout << std::endl;
                }
		#pragma omp barrier
		#pragma omp single
                {
                    gpuGetVelocities();
                }
		#pragma omp barrier
		#pragma omp single
                {
                    //ReSort(octree.root, octree.rootInfo);
                    /*if(ReSort(octree.root, octree.rootInfo))
                     *       RebuildTree(number_particles, precision, maximum_extent, minimum_extent);
                     * */
                }
            }
//             std::cout << "gpu_velocities = ";std::copy(m_gpu_velocity.begin(), m_gpu_velocity.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << std::endl;
	    logger.startTimer("copyVelocities");
            copyVelocities(v);
	    logger.stopTimer("copyVelocities");
//             std::cout << "hv_velocities = ";std::copy(v, v + 3*m_num_particles, std::ostream_iterator<value_type>(std::cout, " ")); std::cout << std::endl;
        }

        void initData(value_type *x, value_type *v, value_type *f)
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
            std::transform(m_gpu_velocity.begin(), m_gpu_velocity.end(), v, v, std::plus<value_type>());
            std::fill(m_gpu_velocity.begin(), m_gpu_velocity.end(), 0.0);
        }
        void setDomain(value_type domain[2][3])
        {
            value_type tmp[3] = {domain[0][0] + domain[1][0], domain[0][1] + domain[1][0], domain[0][2] + domain[1][0]};
            m_extents[0] = 0, m_extents[1] = 0;
            for (int i = 0; i < 3; ++i)
            {
                if (m_extents[0] < tmp[i])
                    m_extents[0] = tmp[i];
                if (m_extents[1] > tmp[i])
                    m_extents[1] = tmp[i];
            }

        }

        void allPairs()
        {
            AllPairs(m_num_particles, 0, m_delta);

            std::cout << "ap_velocities = ";std::copy(m_gpu_velocity.begin(), m_gpu_velocity.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << std::endl;
        }
        void setDelta(value_type delta) { m_delta = delta; }
};


#endif


