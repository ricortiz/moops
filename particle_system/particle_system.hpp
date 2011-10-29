#ifndef PARTICLE_SYSTEM_HPP
#define PARTICLE_SYSTEM_HPP

#include<limits>
#include<algorithm>
#include "particle_system/particle_system_storage.hpp"
#include "particle_system/vtk_particle_system_storage.hpp"
#include "particle_system/particle.hpp"

template<typename _value_type,
         typename _sdc_type,
         size_t _num_particles,
         int immersed_structure_type = PSYS::SURFACE,
         typename storage_type = vtk_particle_system_storage<_value_type,Particle<_value_type>,immersed_structure_type>
         >
class ParticleSystem
{
    public:
        typedef Particle<_value_type> particle_type;
        typedef _value_type          value_type;
        typedef _sdc_type            particle_integrator_type;

        enum
        {
            num_particles = _num_particles
        };
    private:
        storage_type          m_storage;
        value_type            m_time;         ///< Current time.

    public:

        ParticleSystem() : m_time(0.0), m_storage(num_particles) {}
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
            for (int i = 0; i < num_particles; i++)
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
                center[i] /= num_particles;
                max[i] = std::abs(max[i]-center[i])+.001;
                min[i] = std::abs(min[i]-center[i])+.001;
                extent[i] = std::max(max[i],min[i]);
            }
        }

        void clear() { m_time = 0.; }
        void clear_forces()
        {
            size_t size = 3*num_particles;
            value_type *f = forces();
            std::fill(f,f+size,0.0);
        }
        void clear_velocities()
        {
            size_t size = 3*num_particles;
            value_type *v = velocities();
            std::fill(v,v+size,0.0);
        }

        std::size_t particles_size() const
        {
            return num_particles;
        }

        storage_type *storage() { return &m_storage; }
};

// template<typename _value_type, typename _particle_integrator_type, size_t _num_particles>
// struct particle_system_traits<ParticleSystem<_value_type,_particle_integrator_type,_num_particles> >
// {
//     enum
//     {
//         num_particles = _num_particles
//     };
//     typedef _value_type value_type;
//     typedef Particle<value_type> particle_type;
//     typedef _particle_integrator_type particle_integrator_type;
// };


#endif
