#ifndef SINE_GEOMETRY_HPP
#define SINE_GEOMETRY_HPP
/****************************************************************************
** MOOPS -- Modular Object Oriented Particle Simulator
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
#include<vector>
#include <cmath>
#include "base_geometry.hpp"

template<typename value_type>
class SineGeometry : public BaseGeometry<SineGeometry<value_type> >
{

    private:
        value_type       m_length;
        value_type       m_speed;
        value_type       m_amplitude;
        value_type       m_pitch;
        value_type       m_tail_radius;
        value_type       m_head_radius;
        size_t           m_dims[4];
        size_t           m_num_particles;
        value_type       m_x0[3];

    public:
        SineGeometry()
            :       m_length(4.0),
                    m_speed(.001),
                    m_amplitude(.25),
                    m_pitch(4.1),
                    m_tail_radius(.05),
                    m_num_particles(0) {}
        SineGeometry(value_type x0[3], size_t M = 6, size_t N = 100, value_type length = 4.0, value_type speed = .001, value_type amplitude = .25, value_type pitch = 4.1, value_type radius = .05)
            :       m_length(length),
                    m_speed(speed),
                    m_amplitude(amplitude),
                    m_pitch(pitch),
                    m_tail_radius(radius),
                    m_num_particles(0)
        {
            m_dims[0] = M;
            m_dims[1] = N;
            m_dims[2] = 12;
            m_dims[3] = 21;
            m_x0[0] = x0[0];
            m_x0[1] = x0[1];
            m_x0[2] = x0[2];
            m_num_particles = m_dims[0] * m_dims[1] + m_dims[2] * (m_dims[3] - 1) + 1;
        }

        void getDimensions(size_t &Mt, size_t &Nt, size_t &Mh, size_t &Nh)
        {
            Mt = m_dims[0];
            Nt = m_dims[1];
            Mh = m_dims[2];
            Nh = m_dims[3];
        }
        void getDimensions(size_t &Mt, size_t &Nt)
        {
            Mt = m_dims[0];
            Nt = m_dims[1];
        }
        void getHeadRadius(value_type &radius) { radius = m_head_radius; }
        void getTailRadius(value_type &radius) { radius = m_tail_radius; }
        void getWaveSpeed(value_type &speed) { speed = m_speed; }
        void getX0(value_type *x0) { x0 = m_x0; }
        void getLength(value_type &length) { length = m_length; }
        void getTailAmplitude(value_type &amplitude) { amplitude = m_amplitude; }
        void getTailPitch(value_type &pitch) { pitch = m_pitch; }
        void setDimensions(size_t dims[4])
        {
            m_dims[0] = dims[0];
            m_dims[1] = dims[1];
            m_dims[2] = dims[2];
            m_dims[3] = dims[3];
            m_num_particles = dims[0] * dims[1] + dims[2] * (dims[3] - 1) + 1;
        }
        void setDimensions(size_t Mt, size_t Nt, size_t Mh, size_t Nh)
        {
            m_dims[0] = Mt;
            m_dims[1] = Nt;
            m_dims[2] = Mh;
            m_dims[3] = Nh;
            m_num_particles = Mt * Nt + Mh * (Nh - 1) + 1;
        }
        void setHeadRadius(value_type radius) { m_head_radius = radius; }
        void setTailRadius(value_type radius) { m_tail_radius = radius; }
        void setWaveSpeed(value_type speed) { m_speed = speed; }
        void setX0(value_type x[3]) { m_x0[0] = x[0]; m_x0[1] = x[1]; m_x0[2] = x[2]; }
        void setX0(value_type x, value_type y, value_type z)
        {
            m_x0[0] = x;
            m_x0[1] = y;
            m_x0[2] = z;
        }
        void setLength(value_type length) { m_length = length; }
        void setTailAmplitude(value_type amplitude) { m_amplitude = amplitude; }
        void setTailPitch(value_type pitch) { m_pitch = pitch; }

        size_t const &numParticles() const { return m_num_particles; }

        /**
         * @brief Initialize the given particle array in to the geometry shape
         *
         * @param particles array of particles
         * @param grid grid map, maps particle to a grid point
         **/
        template < typename particle_type>
        void init(particle_type *particles)
        {
            value_type dtheta = 2 * M_PI / m_dims[2];
            value_type dalpha = 1.0 / m_dims[3];
            size_t idx = 1;
            particles[0].position[0] = particles[0].position[1] = particles[0].position[2] = 0;
            particles[0].i = 0;
            particles[0].j = 0;
            for(size_t i = 1; i < m_dims[3]; ++i)
                for(size_t j = 0; j < m_dims[2]; ++j, ++idx)
                {
                    headPoint(i, j,particles[idx], dtheta, dalpha);
                    particles[idx].i = i;
                    particles[idx].j = j;
                }
            value_type tail_translation_factor = m_length / 4.0 - .5 * dalpha;
            dtheta = 2 * M_PI / m_dims[0];
            dalpha = 2 * M_PI / m_dims[1];
            for(size_t i = 0; i < m_dims[1]; ++i)
                for(size_t j = 0; j < m_dims[0]; ++j, ++idx)
                {
                    surfacePoint(i, j, 0.0, particles[idx], dtheta, dalpha);
                    particles[idx].position[1] += tail_translation_factor;
                    particles[idx].i = i;
                    particles[idx].j = j;
                }
        }

        inline void surfacePoint(size_t i, size_t j, value_type t, value_type *position, value_type dtheta, value_type dalpha)
        {
            value_type s = i * dalpha;
            value_type x[2] = {0}, normal[2] = {0};
            getTailFrame(s, t, x, normal);
            value_type theta = j * dtheta;
            position[0] = m_tail_radius * normal[0] * std::cos(theta) + x[0];
            position[1] = m_tail_radius * normal[1] * std::cos(theta) + x[1];
            position[2] = m_tail_radius * std::sin(theta)                   ;
        }

        template<typename particle_type>
        void surfacePoint(size_t i, size_t j, value_type t, particle_type &particle, value_type dtheta, value_type dalpha)
        {
            surfacePoint(i, j, t, particle.position, dtheta, dalpha);
        }

        template<typename particle_type>
        void headPoint(size_t i, size_t j,particle_type &particle, value_type dtheta, value_type dalpha)
        {
            value_type R = i < 7 ? m_head_radius * (std::cos((7 - i) * .5 * M_PI / 7)) : i > 13 ? m_head_radius * (std::cos((i - 13) * .5 * M_PI / 6)) : m_head_radius;
            if(i == m_dims[3] - 2 || i == m_dims[3] - 1)
                R = m_tail_radius;
            value_type s = (i - .8) * dalpha;
            /*centerline and tangent*/

            value_type normal = 1.0;
            value_type binormal = 1.0;
            value_type theta = j * dtheta;
            particle.position[0] = R * normal * std::cos(theta);
            particle.position[1] = s;
            particle.position[2] = R * binormal * std::sin(theta);
        }

        template<typename array_type>
        void getTailConnections(array_type &col_ptr, array_type &col_idx, int offset = 0)
        {
            int head_offset = m_dims[2] * (m_dims[3] - 1) + 1;
            for(size_t j = 0; j < m_dims[1]; ++j)
                for(size_t i = 0; i < m_dims[0]; ++i)
                {
                    this->add_plane_connections(i, j, m_dims[0], m_dims[1], col_idx, head_offset + offset);
                    this->add_cylinder_connections(i, j, m_dims[0], m_dims[1], col_idx, head_offset + offset);
                    reinforce_tail(i, j, m_dims[0], m_dims[1], col_idx, head_offset + offset);
                    col_ptr.push_back(col_idx.size());
                }
        }

        template<typename array_type>
        void getHeadConnections(array_type &col_ptr, array_type &col_idx, int offset = 0)
        {
            for(size_t i = 0; i < m_dims[2]; ++i)
                col_idx.push_back(i + 1 + offset);
            col_ptr.push_back(col_idx.size());
            for(size_t j = 0; j < m_dims[3] - 1; ++j)
                for(size_t i = 0; i < m_dims[2]; ++i)
                {
                    this->add_plane_connections(i, j, m_dims[2], m_dims[3] - 1, col_idx, 1 + offset);
                    this->add_cylinder_connections(i, j, m_dims[2], m_dims[3] - 1, col_idx, 1 + offset);
                    if(j == m_dims[3] - 2)
                    {
                        size_t head_offset = m_dims[2] * (m_dims[3] - 1) + 1;
                        for(size_t k = head_offset; k < head_offset + m_dims[0]; ++k)
                            col_idx.push_back(k + offset);
                    }
                    col_ptr.push_back(col_idx.size());
                }
        }

        template<typename array_type>
        void getConnections(array_type &col_ptr, array_type &col_idx, int offset = 0)
        {
            getHeadConnections(col_ptr,col_idx,offset);
            getTailConnections(col_ptr,col_idx,offset);
        }

        template<typename array_type>
        inline void reinforce_tail(int i, int j, int M, int N, array_type &col_idx, size_t offset = 0)
        {
            int connections[1][2] = {{i+M/2,j}};
            col_idx.push_back(connections[0][1]*M + (connections[0][0] + M) % M + offset);
        }

        template <typename particle_type>
        void applyRotation(particle_type *particles, value_type Q[3][3])
        {
            size_t i;
            #pragma omp parallel for private(i) shared(particles,Q)
            for(i = 0; i < m_num_particles; ++i)
            {
                particles[i].position[0] = Q[0][0]*particles[i].position[0]+Q[0][1]*particles[i].position[1]+Q[0][2]*particles[i].position[2];
                particles[i].position[1] = Q[1][0]*particles[i].position[0]+Q[1][1]*particles[i].position[1]+Q[1][2]*particles[i].position[2];
                particles[i].position[2] = Q[2][0]*particles[i].position[0]+Q[2][1]*particles[i].position[1]+Q[2][2]*particles[i].position[2];
            }
        }


    private:
        value_type amplitude_cubic(const value_type &s)
        {
            value_type value = -m_amplitude * 0.03125 * (-70 + 84 * s - 35 * s*s + 5 * s*s*s) * s*s*s*s;
            return s >= 2 ? m_amplitude : s <= 0 ? 0. : value;
        }

        value_type damplitude_cubic(const value_type &s)
        {
            value_type value = -m_amplitude * 0.03125 * 35 * s*s*s * (s - 2)*(s - 2)*(s - 2);
            return (s >= 2 || s <= 0) ? 0. : value;
        }

        void getTailFrame(value_type s, value_type t, value_type x[2], value_type normal[2])
        {
            value_type A = amplitude_cubic(s);
            value_type DA = damplitude_cubic(s);
            value_type Dx[2];

            x[0] = A * std::cos(m_pitch * s - t * m_speed);
            x[1] = s;

            /*tangential vector: derivatives*/
            Dx[0] = DA * std::cos(m_pitch * s - t * m_speed) - m_pitch * A * std::sin(m_pitch * s - t * m_speed);
            Dx[1] = 1;

            value_type norm = std::sqrt(Dx[0] * Dx[0] + Dx[1] * Dx[1]);
            Dx[0] /= norm;
            Dx[1] /= norm;

            /*normal vector*/
            normal[0] = Dx[1];
            normal[1] = -Dx[0];
        }
};

template<typename _value_type>
struct geometry_traits<SineGeometry<_value_type> >
{
    typedef _value_type          value_type;
};

#endif


