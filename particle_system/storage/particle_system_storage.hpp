#ifndef PARTICLE_SYSTEM_STORAGE_HPP
#define PARTICLE_SYSTEM_STORAGE_HPP
/****************************************************************************
 * * MOOPS -- Modular Object Oriented Particle Simulator
 ** Copyright (C) 2011-2012  Ricardo Ortiz <ortiz@unc.edu>
 **
 ** This program is free software: you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation, either version 3 of the License, or
 ** (at your option) any later version.
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ** GNU General Public License for more details.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/
/// @name ParticleSystemStorage - Provides the storage for the particle system.
///
/// @section Description This particle system provides allocation and handling of 
/// 			the simulation data.  The purpose of this class is to provide 
/// 			storage model that is easy to swap when needing.
/// 			

enum storage_shape { SURFACE, VOLUME };

template<typename T, typename particle_type>
struct particle_system_arrays
{
    T positions;
    T velocities;
    T forces;
    particle_type particles;
};

template<typename T, typename particle_type, storage_shape immerse_structure_type>
class ParticleSystemStorage;


template<typename T, typename particle_type>
class ParticleSystemStorage<T,particle_type,SURFACE>
{
    private:
        particle_system_arrays<T*,particle_type*> m_data;
        size_t m_data_size;
        
    public:
        inline explicit ParticleSystemStorage(size_t num_particles)
        {
            m_data_size = num_particles*3;
            m_data.positions = new T[m_data_size];
            m_data.velocities = new T[m_data_size];
            m_data.forces = new T[m_data_size];
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

        size_t data_size() { return m_data_size; }
};

template<typename T, typename particle_type>
class ParticleSystemStorage<T,particle_type,VOLUME>
{
    private:
        particle_system_arrays<T*,particle_type*> m_data;
        size_t m_data_size;
        
    public:
        inline explicit ParticleSystemStorage(size_t num_particles)
        {
            m_data_size = num_particles*3;
            m_data.positions = new T[m_data_size];
            m_data.velocities = new T[m_data_size];
            m_data.forces = NULL;
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
        
        size_t data_size() { return m_data_size; }
};

#endif

