#ifndef PARTICLE_SYSTEM_HPP
#define PARTICLE_SYSTEM_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    ParticleSystem
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @NAME ParticleSystem - Main kernel for the particle system
/// @SECTION Description The ParticleSystem class is the main class for the simulator.  
///			 It owns the storage and provides methods for accessing the data.

#include<limits>
#include<algorithm>

template<typename Derived>
class ParticleSystem
{
    public:
        typedef typename Traits<Derived>::value_type value_type;
        typedef typename Traits<Derived>::particle_type particle_type;
        typedef typename Traits<Derived>::storage_type storage_type;
	

    private:
        value_type            m_time;         ///< Current time.
        storage_type          m_storage;
	size_t 		      m_num_particles;
        value_type            m_domain[2][3]
;
    public:

        ParticleSystem(size_t num_particles) : m_time(0.0), m_storage(num_particles), m_num_particles(num_particles)
        {
        }
        ~ParticleSystem() { };

    public:

        inline value_type & time() { return m_time; }
        inline value_type const & time() const { return m_time; }
        inline value_type *positions()  { return m_storage.positions(); }
        inline value_type const *positions()  const { return m_storage.positions(); }
        inline value_type *velocities() { return m_storage.velocities(); }
        inline value_type const *velocities() const { return m_storage.velocities(); }
        inline value_type *forces()     { return m_storage.forces(); }
        inline value_type const *forces() const { return m_storage.forces(); }
        inline particle_type const *particles() const { return m_storage.particles(); }
        inline particle_type *particles()             { return m_storage.particles(); }


        inline void setDomain()
        {
            particle_type *p = m_storage.particles();
            value_type min[3] = {p[0].position[0], p[0].position[1], p[0].position[2]};
            value_type max[3] = {p[0].position[0], p[0].position[1], p[0].position[2]};
            // center
            m_domain[0][0] = m_domain[0][1] = m_domain[0][2] = 0;
            // extent
            m_domain[1][0] = m_domain[1][1] = m_domain[1][2] = 0;
            for(int i = 0; i < m_num_particles; i++)
            {
                for(int k = 0 ; k < 3; ++k)
                {
                    if(p[i].position[k] > max[k])
                        max[k] = p[i].position[k];
                    if(p[i].position[k] < min[k])
                        min[k] = p[i].position[k];
                    m_domain[0][k] += p[i].position[k];
                }
            }
            
            value_type R0 = 0;
            for(int i = 0; i < 3; ++i)
            {
                m_domain[0][i] /= m_num_particles;
                m_domain[0][i] = int(m_domain[0][i] + .5); // shift center to nearest integer
                R0 = std::max(max[i] - m_domain[0][i], R0);
                R0 = std::max(m_domain[0][i] - min[i], R0);
                max[i] = std::abs(max[i] - m_domain[0][i]) + .001;
                min[i] = std::abs(min[i] - m_domain[0][i]) + .001;
            }
            
            m_domain[1][0] = m_domain[1][1] = m_domain[1][2] = R0 * 1.000001;
        }

        void clear() { m_time = 0.; }
        inline void clear_forces() { std::fill(forces(),forces()+data_size(),0.0); }
        inline void clear_velocities() { std::fill(velocities(),velocities()+data_size(),0.0); }
        inline std::size_t particles_size() { return m_num_particles; }      
        inline std::size_t data_size() { return m_storage.data_size(); }
        storage_type &storage() { return m_storage; }
};

#endif
