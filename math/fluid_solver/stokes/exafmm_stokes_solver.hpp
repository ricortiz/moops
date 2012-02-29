#ifndef EXAFMM_STOKES_SOLVER_HPP
#define EXAFMM_STOKES_SOLVER_HPP

#include <vector>

struct JBody
{
    int         IBODY;                                            //!< Initial body numbering for sorting back
    int         IPROC;                                            //!< Initial process numbering for partitioning back
    unsigned    ICELL;                                            //!< Cell index
    vect        X;                                                //!< Position
    vect        FORCE;                                            //!< Force
};

struct Body : public JBody
{
    vec<3, real> TRG;                                             //!< velocity values
    bool operator<(const Body &rhs) const                         //!< Overload operator for comparing body index
    {
        return this->IBODY < rhs.IBODY;                             //!< Comparison function for body index
    }
};

typedef Body exafmm_jparticle;

#include "fmm/exafmm/include/serialfmm.h"

template < typename value_type>
class ExaFmmStokesSolver
{
protected:
    typedef std::vector<exafmm_particle>        particle_array_type;
    typedef std::vector<Cell>                   box_array_type;
    typedef particle_array_type::iterator       particle_iterator;

    private:
        size_t                  m_num_particles;
        particle_array_type     m_particles;
        box_array_type          m_target_boxes;
        box_array_type          m_source_boxes;
        SerialFMM<Stokes>       m_fmm;

    public:
        ExaFmmStokesSolver(size_t num_particles)
                :
                m_num_particles(num_particles),
                m_particles(num_particles)                
        {
            m_fmm.initialize();            
        }

        ~ExaFmmStokesSolver()
        {
            m_fmm.finalize();
        }
        
        inline void operator()(value_type, value_type *x, value_type *v, value_type *f)
        {
            size_t idx = 0;
            for( particle_iterator p = m_particles.begin(), end = m_particles.end(); p!= end; ++p )
            {
                p->X[0] = x[idx];
                p->X[1] = x[idx+1];
                p->X[2] = x[idx+2];
                
                p->FORCE[0] = f[idx];
                p->FORCE[1] = f[idx+1];
                p->FORCE[2] = f[idx+2];
                p->IBODY = p-m_particles.begin();
                p->IPROC = MPIRANK;
                idx += 3;
            }
            m_fmm.setDomain(m_particles);
            m_source_boxes.clear();
            m_fmm.bottomup(m_particles, m_target_boxes);
            m_source_boxes = m_target_boxes;
            m_fmm.downward(m_target_boxes, m_source_boxes);
            for( particle_iterator p = m_particles.begin(), end = m_particles.end(); p!= end; ++p )
            {
                idx = 3*p->IBODY;
                v[idx] = p->TRG[0];
                v[idx+1] = p->TRG[1];
                v[idx+2] = p->TRG[2];
            }
        }

        void setDelta(value_type delta) { m_delta = delta; }

        void allPairs()
        {
            std::vector<value_type> velocity(3*m_num_particles, 0.0);
            for (size_t i = 0; i < m_num_particles; ++i)
            {
                unsigned i_idx = 3 * octree.particle_idx[i];
                for (size_t j = 0; j < m_num_particles; ++j)
                    computeStokeslet(octree.bodies[i].position, &velocity[i_idx], octree.bodies[j].position, octree.bodies[j].force, m_delta);
            }
            std::cout << "ap_velocities = [";std::copy(velocity.begin(), velocity.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "]" << std::endl;
        }

};


#endif


