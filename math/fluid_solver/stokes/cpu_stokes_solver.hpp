#ifndef CPU_STOKES_SOLVER_HPP
#define CPU_STOKES_SOLVER_HPP

#include <vector>
#include <map>
#include "nbody_cpu/cpu_compute_velocity.hpp"

template<typename value_type>
class CpuStokesSolver
{
    protected:
        typedef std::vector<std::vector<size_t> > map_type;

    private:

        value_type      m_delta;
        map_type        m_map_explicit;
        map_type        m_map_implicit;
        size_t          m_num_particles;

    public:
        CpuStokesSolver(size_t num_particles) : m_num_particles(num_particles) {}
        
        inline void operator() ( value_type t, const value_type *x, value_type *v, const  value_type *f)
        {
            operator() ( t, x, v, x, f, m_num_particles, m_num_particles );
        }

        inline void operator() ( value_type t, const value_type *x, value_type *v, const value_type *y, const value_type *f, size_t num_sources, size_t num_targets )
        {
            size_t size_targets = 3 * num_targets;
            size_t size_sources = 3 * num_sources;
            std::fill(v,v+size_targets,0.0);
            #pragma omp parallel for shared(y,f)
            for ( size_t i = 0; i < size_targets; i += 3 )
                for ( size_t j = 0; j < size_sources; j += 3 )
                    compute_velocity ( &x[i], &v[i], &y[j], &f[j], m_delta );
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
            std::fill(v,v+m_num_particles*3,0.0);
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
            std::fill(v,v+m_num_particles*3,0.0);
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

