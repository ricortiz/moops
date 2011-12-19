#ifndef PARTICLE_SYSTEM_STORAGE_HPP
#define PARTICLE_SYSTEM_STORAGE_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    ParticleSystemStorage
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name ParticleSystemStorage - Provides the storage for the particle system.
///
/// @section Description This particle system provides allocation and handling of 
/// 			the simulation data.  The purpose of this class is to provide 
/// 			storage model that is easy to swap when needing.
/// 			

namespace PSYS
{
    enum { SURFACE, VOLUME };
}


template<typename T, typename particle_type>
struct particle_system_arrays
{
    T *positions;
    T *velocities;
    T *forces;
    particle_type *particles;
};

template<typename T, typename particle_type, int immerse_structure_type>
class ParticleSystemStorage;


template<typename T, typename particle_type>
class ParticleSystemStorage<T,particle_type,PSYS::SURFACE>
{
    private:
        particle_system_arrays<T,particle_type> m_data;

    public:
        inline explicit ParticleSystemStorage(size_t data_size)
        {
            size_t num_particles = data_size/3;
            m_data.positions = new T[data_size];
            m_data.velocities = new T[data_size];
            m_data.forces = new T[data_size];
            m_data.particles = new particle_type[num_particles];
            for (size_t i = 0, idx = 0; i < num_particles; ++i, idx+=3)
            {
                m_data.particles[i].position = &m_data.positions[idx];
                m_data.particles[i].velocity = &m_data.velocities[idx];
                m_data.particles[i].force = &m_data.forces[idx];
            }
        }

        ~ParticleSystemStorage()
        {
            delete [] m_data.positions;
            delete [] m_data.velocities;
            delete [] m_data.forces;
            delete [] m_data.particles;
        }
        inline void swap(ParticleSystemStorage& other) { std::swap(m_data,other.m_data); }

        inline const T *position(size_t i) const { return m_data.particle[i].position; }
        inline T *position(size_t i)             { return m_data.particle[i].position; }

        inline const T *velocity(size_t i) const { return m_data.particle[i].velocity; }
        inline T *velocity(size_t i)             { return m_data.particle[i].velocity; }

        inline const T *force(size_t i) const { return m_data.particle[i].force; }
        inline T *force(size_t i)             { return m_data.particle[i].force; }

        inline const particle_type *particles(size_t i) const { return m_data.particles[i]; }
        inline particle_type *particles(size_t i)             { return m_data.particles[i]; }

        inline const T *positions() const { return m_data.positions; }
        inline T *positions()             { return m_data.positions; }

        inline const T *velocities() const { return m_data.velocities; }
        inline T *velocities()             { return m_data.velocities; }

        inline const T *forces() const { return m_data.forces; }
        inline T *forces()             { return m_data.forces; }

        inline const particle_type *particles() const { return m_data.particles; }
        inline particle_type *particles()             { return m_data.particles; }
};

template<typename T, typename particle_type>
class ParticleSystemStorage<T,particle_type,PSYS::VOLUME>
{
    private:
        particle_system_arrays<T,particle_type> m_data;

    public:
        inline explicit ParticleSystemStorage(size_t data_size)
        {
            size_t num_particles = data_size/3;
            m_data.positions = new T[data_size];
            m_data.velocities = new T[data_size];
            m_data.forces = 0;
            m_data.particles = new particle_type[num_particles];
            for (size_t i = 0, idx = 0; i < num_particles; ++i, idx+=3)
            {
                m_data.particles[i].position = &m_data.positions[idx];
                m_data.particles[i].velocity = &m_data.velocities[idx];
            }
        }

        ~ParticleSystemStorage()
        {
            delete [] m_data.positions;
            delete [] m_data.velocities;
            delete [] m_data.particles;
        }
        inline void swap(ParticleSystemStorage& other) { std::swap(m_data,other.m_data); }

        inline const T *position(size_t i) const { return m_data.particle[i].position; }
        inline T *position(size_t i)             { return m_data.particle[i].position; }

        inline const T *velocity(size_t i) const { return m_data.particle[i].velocity; }
        inline T *velocity(size_t i)             { return m_data.particle[i].velocity; }

        inline const particle_type *particles(size_t i) const { return m_data.particles[i]; }
        inline particle_type *particles(size_t i)             { return m_data.particles[i]; }
        
        inline const T *force(size_t i) const { return m_data.particle[i].force; }
        inline T *force(size_t i)             { return m_data.particle[i].force; }

        inline const T *positions() const { return m_data.positions; }
        inline T *positions()             { return m_data.positions; }
        
        inline const T *forces() const { return m_data.forces; }
        inline T *forces()             { return m_data.forces; }

        inline const T *velocities() const { return m_data.velocities; }
        inline T *velocities()             { return m_data.velocities; }

        inline const particle_type *particles() const { return m_data.particles; }
        inline particle_type *particles()             { return m_data.particles; }
};


#endif

