#ifndef HYBRID_FMM_STOKES_SOLVER_HPP
#define HYBRID_FMM_STOKES_SOLVER_HPP

#include <vector>
#include <map>
extern "C" {
#include "math/fluid_solver/stokes/fmm/hybrid/hybridfmm.h"
}
template<typename value_type>
class HybridFmmStokesSolver
{
    protected:
        typedef std::vector<std::vector<size_t> > map_type;

    private:
        value_type m_delta;
        map_type m_map_explicit;
        map_type m_map_implicit;
        size_t m_num_particles;

    public:
        HybridFmmStokesSolver(size_t num_particles) : m_num_particles(num_particles)
        {
            octree.maxBodiesPerNode = 50;
            octree.numParticles = num_particles;
            
            int precision = 6;

            
            
        }
        
        inline void operator() ( value_type t, const value_type *x, value_type *v, const  value_type *f)
        {
            operator() ( t, x, v, x, f, m_num_particles, m_num_particles );
        }

        inline void operator() ( value_type t, const value_type *x, value_type *v, const value_type *y, const value_type *f, size_t num_sources, size_t num_targets )
        {
            size_t size_targets = 3 * num_targets;
            size_t size_sources = 3 * num_sources;
            
        }

        inline void Implicit( value_type t, const value_type *x, value_type *v, const value_type *f )
        {
            Implicit( t, x, v, x, f );
        }

        inline void Explicit(value_type t, const value_type *x, value_type *v, const value_type *f )
        {
            Explicit  ( t, x, v, x, f );
        }

        inline void Implicit ( value_type t, const value_type *x, value_type *v, const value_type *y, const value_type *f )
        {
            typedef std::vector<size_t>::iterator iterator;
            #pragma omp parallel for shared(y,f)
            for ( size_t i = 0; i < m_num_particles; ++i )
                for(iterator j = m_map_implicit[i].begin(), end = m_map_implicit[i].end(); j != end; ++j)
                {
                    size_t yidx = *j;
                    compute_velocity ( &x[3 * i], &v[3 * i], &y[yidx], &f[yidx], m_delta );
                }
        }

        inline void Explicit ( value_type t, const value_type *x, value_type *v, const value_type *y, const value_type *f )
        {
            typedef std::vector<size_t>::iterator iterator;
            #pragma omp parallel for shared(y,f)
            for ( size_t i = 0; i < m_num_particles; ++i)
                for(iterator j = m_map_explicit[i].begin(), end = m_map_explicit[i].end(); j != end; ++j)
                {
                    size_t yidx = *j;
                    compute_velocity ( &x[3 * i], &v[3 * i], &y[yidx], &f[yidx], m_delta );
                }
        }

        template<typename spring_system_type>
        inline void initMaps ( spring_system_type &spring_system )
        {
            typedef typename spring_system_type::particle_type   particle_type;
            typedef typename spring_system_type::spring_lut_type spring_map_type;
            m_map_explicit.resize(m_num_particles);
            m_map_implicit.resize(m_num_particles);
            particle_type *particles = spring_system.derived().particles();
            for ( size_t i = 0; i < m_num_particles; ++i)
                for ( size_t j = 0, yidx = 0; j < m_num_particles; ++j, yidx += 3 )
                    if ( spring_system.existSpring ( &particles[i], &particles[j] ) )
                        m_map_implicit[i].push_back(yidx);
                    else
                        m_map_explicit[i].push_back(yidx);
        }
        
        void setDelta(value_type delta) { m_delta = delta; }
};


#endif

