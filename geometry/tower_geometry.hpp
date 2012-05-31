#ifndef TOWER_GEOMETRY_HPP
#define TOWER_GEOMETRY_HPP
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

#include "base_geometry.hpp"

template<typename value_type>
class TowerGeometry : public BaseGeometry<TowerGeometry<value_type> >
{
    
private:
    value_type       m_length;
    value_type       m_radius;
    size_t           m_dims[2];
    size_t           m_num_particles;
    value_type       m_x0[3];
    
public:
    TowerGeometry() : m_length(20.0), m_num_particles(0), m_radius(1.0)
    {
        m_x0[0] = m_x0[1] = m_x0[2] = 0.0;
        m_dims[0] = m_dims[1] = 0;
    }
    TowerGeometry(value_type x0[3], size_t M = 6, size_t N = 20, value_type length = 20.0, value_type radius = .05)
        : m_length(length), m_num_particles(0), m_radius(radius)
    {
        m_dims[0] = M;
        m_dims[1] = N;
        m_x0[0] = x0[0];
        m_x0[1] = x0[1];
        m_x0[2] = x0[2];
        m_num_particles = m_dims[0] * m_dims[1];
    }

    void getDimensions(size_t &M, size_t &N)
    {
        M = m_dims[0];
        N = m_dims[1];
    }
    void getRadius(value_type &radius)  { radius = m_radius; }
    void getX0(value_type *x0)          { x0 = m_x0; }
    void getLength(value_type &length)  { length = m_length; }
    void setDimensions(size_t dims[2])
    {
        m_dims[0] = dims[0];
        m_dims[1] = dims[1];
        m_num_particles = dims[0] * dims[1];
    }
    void setDimensions(size_t M, size_t N)
    {
        m_dims[0] = M;
        m_dims[1] = N;
        m_num_particles = M * N;
    }
    void setRadius(value_type radius) { m_radius = radius; }
    void setX0(value_type x[3]) { m_x0[0] = x[0]; m_x0[1] = x[1]; m_x0[2] = x[2]; }
    void setX0(value_type x, value_type y, value_type z)
    {
        m_x0[0] = x;
        m_x0[1] = y;
        m_x0[2] = z;
    }
    void setLength(value_type length) { m_length = length; }
    
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
        value_type dtheta = 2 * M_PI / m_dims[0];
        value_type ds = m_length/m_dims[1];
        size_t idx = 1;
        for(size_t i = 0; i < m_dims[1]; ++i)
            for(size_t j = 0; j < m_dims[0]; ++j, ++idx)
            {
                surfacePoint(i, j, particles[idx], dtheta, ds);
                particles[idx].i = i;
                particles[idx].j = j;
            }
    }
    
    inline void surfacePoint(size_t i, size_t j, value_type *position, value_type dtheta, value_type ds)
    {
        value_type x[2] = {0}, normal[2] = {0};
        value_type theta = j * dtheta;
        position[0] = m_radius * std::cos(theta) + x[0];
        position[1] = m_radius * std::sin(theta) + x[1];
        position[2] = i * ds + x[2];
    }
    
    template<typename particle_type>
    void surfacePoint(size_t i, size_t j, particle_type &particle, value_type dtheta, value_type ds)
    {
        surfacePoint(i, j, particle.position, dtheta, ds);
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
    
};

template<typename _value_type>
struct geometry_traits<TowerGeometry<_value_type> >
{
    typedef _value_type          value_type;
};

#endif


