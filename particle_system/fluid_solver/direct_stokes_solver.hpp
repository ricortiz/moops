#ifndef DIRECT_STOKES_SOLVER_HPP
#define DIRECT_STOKES_SOLVER_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    DirectStokesSolver
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name DirectStokesSolver - Solves the Stokes problem using a direct solver.
/// @section Description DirectStokesSolver solves an n-body problem and updates the
///                      velocities of the particle system.
/// @section See also

#include <map>
#include "math/fluid_solver/stokes/nbody_cpu/cpu_compute_velocity.hpp"

template<typename particle_system_type>
class DirectStokesSolver
{
    protected:
        typedef typename particle_system_type::value_type value_type;
        typedef typename particle_system_type::particle_type particle_type;
        typedef std::map<size_t,std::vector<size_t> > map_type;

    private:

        value_type m_delta;
        size_t m_num_particles;
        map_type m_map_explicit;
        map_type m_map_implicit;

    public:
        DirectStokesSolver() : m_num_particles ( 0 ) {}
        DirectStokesSolver ( size_t num_particles ) : m_num_particles ( num_particles ) {}

        inline void operator() ( const value_type *x, value_type *v,const  value_type *f )
        {
            operator() ( x, v, x, f, m_num_particles, m_num_particles );
        }
        
        inline void operator() ( const value_type *x, value_type *v, const value_type *y, const value_type *f, size_t num_sources, size_t num_targets )
        {
            size_t size_targets = 3 * num_targets;
            size_t size_sources = 3 * num_sources;
            for ( size_t i = 0; i < size_targets; i += 3 )
                for ( size_t j = 0; j < size_sources; j += 3 )
                    compute_velocity ( &x[i], &v[i], &y[j], &f[j], m_delta );
        }
        
        inline void ImplicitOperator( const value_type *x, value_type *v, const value_type *f )
        {
            ImplicitOperator( x, v, x, f );
        }
        
        inline void ExplicitOperator(const  value_type *x, value_type *v, const value_type *f )
        {
            ExplicitOperator  ( x, v, x, f );
        }
        
        inline void ImplicitOperator ( const value_type *x, value_type *v, const value_type *y, const value_type *f )
        {
            typedef std::vector<size_t>::iterator iterator;
            for ( size_t i = 0, xidx = 0; i < m_num_particles; ++i, xidx += 3 )
                for(iterator j = m_map_implicit[i].begin(), end = m_map_implicit[i].end(); j != end; ++j)
                {
                    size_t yidx = *j;
                    compute_velocity ( &x[xidx], &v[xidx], &y[yidx], &f[yidx], m_delta );
                }
        }
        
        inline void ExplicitOperator ( const value_type *x, value_type *v, const value_type *y, const value_type *f )
        {
            typedef std::vector<size_t>::iterator iterator;            
            for ( size_t i = 0, xidx = 0; i < m_num_particles; ++i, xidx += 3 )
                for(iterator j = m_map_explicit[i].begin(), end = m_map_explicit[i].end(); j != end; ++j)
                {
                    size_t yidx = *j;
                    compute_velocity ( &x[xidx], &v[xidx], &y[yidx], &f[yidx], m_delta );
                }
        }

        template<typename spring_system_type>
        inline void initMaps ( spring_system_type &spring_system )
        {
            typedef typename spring_system_type::spring_lut_type spring_map_type;
            particle_type *particles = spring_system.particles();
            
            for ( size_t i = 0, xidx = 0; i < m_num_particles; ++i, xidx += 3)
                for ( size_t j = 0, yidx = 0; j < m_num_particles; ++j, yidx += 3 )
                    if ( spring_system.exist_spring ( &particles[i], &particles[j] ) )
                        m_map_implicit[i].push_back(yidx);
                    else
                        m_map_explicit[i].push_back(yidx);
        }
        
        
        value_type const &delta() const { return m_delta; }
        value_type &delta() { return m_delta; }
};


#endif

