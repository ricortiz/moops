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

/** \class ParticleSystemStorage
 *  \ingroup ParticleSystem_Module
 *  \brief This particle system provides allocation and handling of the simulation data.
 *
 * The purpose of this class is to provide storage model that is easy to swap when needed.
 *
 * \tparam T is the data type
 * \tparam P is the particle wrapper type
 * \tparam immerse_structure_type is the type of immersed particles, free moving or constrained(by spring, for example)
 */
enum storage_shape { SURFACE, VOLUME };

template<typename T, typename P>
struct particle_system_arrays
{
    T positions;
    T velocities;
    T forces;
    P particles;
};

template<typename T, typename P, storage_shape immerse_structure_type>
class ParticleSystemStorage;


template<typename T, typename P>
class ParticleSystemStorage<T,P,SURFACE>
{
    public:
        typedef std::vector<T> data_vector_type;
        typedef typename data_vector_type::iterator data_iterator;
        typedef std::vector<P> particle_vector_type;
        typedef typename particle_vector_type::iterator particle_iterator;

    private:
        particle_system_arrays<data_vector_type,particle_vector_type> m_data;
        size_t m_data_size;

    public:
        inline explicit ParticleSystemStorage(size_t num_particles)
        {
            m_data_size = num_particles*3;
            m_data.positions.resize(m_data_size);
            m_data.velocities.resize(m_data_size);
            m_data.forces.resize(m_data_size);
            m_data.particles.resize(num_particles);
            for(size_t i = 0, idx = 0; i < num_particles; ++i, idx+=3)
            {
                m_data.particles[i].position = &m_data.positions[idx];
                m_data.particles[i].velocity = &m_data.velocities[idx];
                m_data.particles[i].force = &m_data.forces[idx];
            }
        }

        ~ParticleSystemStorage()
        {
        }
        inline void swap(ParticleSystemStorage& other) { std::swap(m_data,other.m_data); }

        inline data_iterator positions_begin() { return m_data.positions.begin(); }
        inline data_iterator positions_end()   { return m_data.positions.end(); }

        inline data_iterator velocities_begin() { return m_data.velocities.begin(); }
        inline data_iterator velocities_end()   { return m_data.velocities.end(); }

        inline data_iterator forces_begin() { return m_data.forces.begin(); }
        inline data_iterator forces_end()   { return m_data.forces.end(); }

        inline particle_iterator particles_begin() { return m_data.particles.begin(); }
        inline particle_iterator particles_end()   { return m_data.particles.end(); }

        inline const T *positions() const { return &m_data.positions[0]; }
        inline T *positions()             { return &m_data.positions[0]; }

        inline const T *velocities() const { return &m_data.velocities[0]; }
        inline T *velocities()             { return &m_data.velocities[0]; }

        inline const T *forces() const { return &m_data.forces[0]; }
        inline T *forces()             { return &m_data.forces[0]; }

        inline const P *particles() const { return &m_data.particles[0]; }
        inline P *particles()             { return &m_data.particles[0]; }

        size_t data_size() { return m_data_size; }

        void clearData()
        {
            std::fill(m_data.positions.begin(),m_data.positions.end(),T(0));
            std::fill(m_data.velocities.begin(),m_data.velocities.end(),T(0));
            std::fill(m_data.forces.begin(),m_data.forces.end(),T(0));
        }
};

template<typename T, typename P>
class ParticleSystemStorage<T,P,VOLUME>
{
    protected:
        typedef std::vector<T> data_vector_type;
        typedef typename data_vector_type::iterator data_iterator;
        typedef std::vector<P> particle_vector_type;
        typedef typename particle_vector_type::iterator particle_iterator;
    private:
        particle_system_arrays<std::vector<T>,std::vector<P> > m_data;
        size_t m_data_size;

    public:
        inline explicit ParticleSystemStorage(size_t num_particles)
        {
            m_data_size = num_particles*3;
            m_data.positions.resize(m_data_size);
            m_data.velocities.resize(m_data_size);
            m_data.particles.resize(num_particles);
            for(size_t i = 0, idx = 0; i < num_particles; ++i, idx+=3)
            {
                m_data.particles[i].position = &m_data.positions[idx];
                m_data.particles[i].velocity = &m_data.velocities[idx];
            }
        }

        ~ParticleSystemStorage()
        {
        }

        inline void swap(ParticleSystemStorage& other) { std::swap(m_data,other.m_data); }

        inline data_iterator positions_begin() { return m_data.positions.begin(); }
        inline data_iterator positions_end()   { return m_data.positions.end(); }

        inline data_iterator velocities_begin() { return m_data.velocities.begin(); }
        inline data_iterator velocities_end()   { return m_data.velocities.end(); }

        inline particle_iterator particles_begin() { return m_data.particles.begin(); }
        inline particle_iterator particles_end()   { return m_data.particles.end(); }

        inline const T *positions() const { return m_data.positions; }
        inline T *positions()             { return m_data.positions; }

        inline const T *forces() const { return m_data.forces; }
        inline T *forces()             { return m_data.forces; }

        inline const T *velocities() const { return m_data.velocities; }
        inline T *velocities()             { return m_data.velocities; }

        inline const P *particles() const { return m_data.particles; }
        inline P *particles()             { return m_data.particles; }

        size_t data_size() { return m_data_size; }

        void clearData()
        {
            std::fill(m_data.positions.begin(),m_data.positions.end(),T(0));
            std::fill(m_data.velocities.begin(),m_data.velocities.end(),T(0));
        }
};

#endif

