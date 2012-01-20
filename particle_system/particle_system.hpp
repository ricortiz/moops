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
#include "particle_system/particle_system_storage.hpp"
#include "particle_system/vtk_particle_system_storage.hpp"
#include "particle_system/particle.hpp"

template<typename _value_type,
         int immersed_structure_type = PSYS::SURFACE,
         typename storage_type = vtkParticleSystemStorage<_value_type,Particle<_value_type>,immersed_structure_type>
         >
class ParticleSystem
{
    public:
        typedef Particle<_value_type> particle_type;
        typedef _value_type          value_type;

    private:
        storage_type          m_storage;
        value_type            m_time;         ///< Current time.
        size_t 		      m_data_size;
	size_t 		      m_num_particles;

    public:

        ParticleSystem(size_t data_size) : m_time(0.0), m_storage(data_size), m_data_size(data_size), m_num_particles(data_size/3)
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


        inline void get_dimensions(value_type center[], value_type extent[])
        {
            value_type highest = std::numeric_limits<value_type>::infinity();
            value_type lowest = std::numeric_limits<value_type>::min();
            value_type min[3] = {highest,highest,highest};
            value_type max[3] = {lowest,lowest,lowest};
            center[0] = center[1] = center[2] = 0;
            extent[0] = extent[1] = extent[2] = 0;
            for (int i = 0; i < m_num_particles; i++)
            {
                value_type *position = m_storage.position(i);
                for (int k = 0 ; k < 3; ++k)
                {
                    if (position[k] > max[k])
                        max[k] = position[k];
                    if (position[k] < min[k])
                        min[k] = position[k];
                    center[k] += position[k];
                }
            }

            for (int i = 0; i < 3; ++i)
            {
                center[i] /= m_num_particles;
                max[i] = std::abs(max[i]-center[i])+.001;
                min[i] = std::abs(min[i]-center[i])+.001;
                extent[i] = std::max(max[i],min[i]);
            }
        }

        void clear() { m_time = 0.; }
        void clear_forces()
        {
            size_t size = 3*m_num_particles;
            value_type *f = forces();
            std::fill(f,f+size,0.0);
        }
        void clear_velocities()
        {
            size_t size = 3*m_num_particles;
            value_type *v = velocities();
            std::fill(v,v+size,0.0);
        }

        std::size_t particles_size() const
        {
            return m_num_particles;
        }
        
        std::size_t data_size() const
        {
            return m_data_size;
        }
        
        storage_type *storage() { return &m_storage; }
};

#endif
