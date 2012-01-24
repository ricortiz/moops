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

        value_type m_delta;
        map_type m_map_explicit;
        map_type m_map_implicit;

    public:

        inline void operator() ( value_type t, const value_type *x, value_type *v, const  value_type *f, size_t num_particles )
        {
            operator() ( x, v, x, f, num_particles, num_particles );
        }

        inline void operator() ( const value_type *x, value_type *v, const value_type *y, const value_type *f, size_t num_sources, size_t num_targets )
        {
            logger.startTimer("computeStokeslets");
            size_t size_targets = 3 * num_targets;
            size_t size_sources = 3 * num_sources;
            #pragma omp parallel for shared(y,f)
            for ( size_t i = 0; i < size_targets; i += 3 )
                for ( size_t j = 0; j < size_sources; j += 3 )
                    compute_velocity ( &x[i], &v[i], &y[j], &f[j], m_delta );
            logger.stopTimer("computeStokeslets");
        }

        inline void Implicit( const value_type *x, value_type *v, const value_type *f, size_t num_particles )
        {
            Implicit( x, v, x, f,num_particles );
        }

        inline void Explicit(const  value_type *x, value_type *v, const value_type *f, size_t num_particles )
        {
            Explicit  ( x, v, x, f,num_particles );
        }

        inline void Implicit ( const value_type *x, value_type *v, const value_type *y, const value_type *f, size_t num_particles )
        {
            typedef std::vector<size_t>::iterator iterator;
            #pragma omp parallel for shared(y,f)
            for ( size_t i = 0; i < num_particles; ++i )
                for(iterator j = m_map_implicit[i].begin(), end = m_map_implicit[i].end(); j != end; ++j)
                {
                    size_t yidx = *j;
                    compute_velocity ( &x[3 * i], &v[3 * i], &y[yidx], &f[yidx], m_delta );
                }
        }

        inline void Explicit ( const value_type *x, value_type *v, const value_type *y, const value_type *f, size_t num_particles )
        {
            typedef std::vector<size_t>::iterator iterator;
            #pragma omp parallel for shared(y,f)
            for ( size_t i = 0; i < num_particles; ++i)
                for(iterator j = m_map_explicit[i].begin(), end = m_map_explicit[i].end(); j != end; ++j)
                {
                    size_t yidx = *j;
                    compute_velocity ( &x[3 * i], &v[3 * i], &y[yidx], &f[yidx], m_delta );
                }
        }

        template<typename spring_system_type>
        inline void initMaps ( spring_system_type &spring_system, size_t num_particles )
        {
            m_map_explicit.resize(num_particles);
            m_map_implicit.resize(num_particles);
            typedef typename spring_system_type::particle_type particle_type;
            typedef typename spring_system_type::spring_lut_type spring_map_type;
            particle_type *particles = spring_system.particles();
            for ( size_t i = 0; i < num_particles; ++i)
                for ( size_t j = 0, yidx = 0; j < num_particles; ++j, yidx += 3 )
                    if ( spring_system.existSpring ( &particles[i], &particles[j] ) )
                        m_map_implicit[i].push_back(yidx);
                    else
                        m_map_explicit[i].push_back(yidx);
        }
        
        void setDelta(value_type delta) { m_delta = delta; }
};


#endif

