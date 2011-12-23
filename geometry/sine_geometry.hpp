#ifndef SINE_GEOMETRY_HPP
#define SINE_GEOMETRY_HPP

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
        value_type       m_radius;
        size_t           m_dims[4];
        size_t           m_size;
        value_type       m_x0[3];

    public:

        SineGeometry()
            : m_length(4.0), m_speed(.001), m_amplitude(.25), m_pitch(4.1), m_radius(.05), m_size(0) {}
        SineGeometry(value_type x0[3], size_t M = 6, size_t N = 100, value_type length = 4.0, value_type speed = .001, value_type amplitude = .25, value_type pitch = 4.1, value_type radius = .05)
            : m_length(length), m_speed(speed), m_amplitude(amplitude), m_pitch(pitch), m_radius(radius), m_size(0)
        {
            m_dims[0] = M;
            m_dims[1] = N;
            m_dims[2] = 12;
            m_dims[3] = 21;
            m_x0[0] = x0[0];
            m_x0[1] = x0[1];
            m_x0[2] = x0[2];
        }

        value_type &length()                    { return m_length; }
        value_type const &length() const        { return m_length; }
        value_type &speed()                     { return m_speed; }
        value_type const &speed() const         { return m_speed; }
        value_type &amplitude()                 { return m_amplitude; }
        value_type const &amplitude() const     { return m_amplitude; }
        value_type &pitch()                     { return m_pitch; }
        value_type const &pitch() const         { return m_pitch; }
        value_type &radius()                    { return m_radius; }
        value_type const &radius() const        { return m_radius; }

        size_t *get_dimensions() { return m_dims; }
        size_t const *get_dimensions() const { return m_dims; }
        size_t size() { return m_size; }
        value_type *getX0() { return m_x0; }

        void setX0(value_type x, value_type y, value_type z)
        {
            m_x0[0] = x;
            m_x0[1] = y;
            m_x0[2] = z;
        }
        /**
         * @brief Initialize the given particle array in to the geometry shape
         *
         * @param particles array of particles
         * @param grid grid map, maps particle to a grid point
         **/
        template<typename particle_type, typename grid_type>
        void init(particle_type *particles, grid_type &grid)
        {
            value_type dtheta = 2 * M_PI / m_dims[2];
            value_type dalpha = m_length / 4.0 / m_dims[3];
            size_t idx = 1;
            particles[0].position[0] = particles[0].position[1] = particles[0].position[2] = 0;
            particles[0].i = 0;
            particles[0].j = 0;
            m_size++;
//             grid[particles] = std::make_pair(0, 0);
            for(size_t i = 1; i < m_dims[3]; ++i)
                for(size_t j = 0; j < m_dims[2]; ++j, ++idx)
                {
                    head_point(i, j, 0.0, particles[idx], dtheta, dalpha);
                    particles[idx].i = i;
                    particles[idx].j = j;
//                     grid[particles+idx] = std::make_pair(i, j);
                    m_size++;
                }
            value_type tail_translation_factor = m_length / 4.0 - .5 * dalpha;
            dtheta = 2 * M_PI / m_dims[0];
            dalpha = 2 * M_PI / m_dims[1];
            for(size_t i = 0; i < m_dims[1]; ++i)
                for(size_t j = 0; j < m_dims[0]; ++j, ++idx)
                {
                    surface_point(i, j, 0.0, particles[idx], dtheta, dalpha);
                    particles[idx].position[1] += tail_translation_factor;
                    particles[idx].i = i;
                    particles[idx].j = j;
                    grid[particles + idx] = std::make_pair(i, j);
                    m_size++;
                }
            for(size_t i = 0; i < idx; ++i)
            {
                particles[i].position[0] += m_x0[0];
                particles[i].position[1] += m_x0[1];
                particles[i].position[2] += m_x0[2];
            }
        }

        inline void surface_point(size_t i, size_t j, value_type t, value_type *position, value_type dtheta, value_type dalpha)
        {
            value_type s = i * dalpha;
            /*centerline and tangent*/
            value_type x[2] = {0}, normal[2] = {0};
            get_tail_frame(s, t, x, normal);
            value_type theta = j * dtheta;
            position[0] = m_radius * normal[0] * std::cos(theta) + x[0];
            position[1] = m_radius * normal[1] * std::cos(theta) + x[1];
            position[2] = m_radius * std::sin(theta)                   ;
        }
        template<typename particle_type>
        void surface_point(size_t i, size_t j, value_type t, particle_type &particle, value_type dtheta, value_type dalpha)
        {
            surface_point(i, j, t, particle.position, dtheta, dalpha);
        }
        template<typename particle_type>
        void head_point(size_t i, size_t j, value_type t, particle_type &particle, value_type dtheta, value_type dalpha)
        {
            value_type scale = 5.0 * m_radius;
            value_type R = i < 7 ? scale * (std::cos((7 - i) * .5 * M_PI / 7)) : i > 13 ? scale * (std::cos((i - 13) * .5 * M_PI / 6)) : scale;
            if(i == m_dims[3] - 2 || i == m_dims[3] - 1)
                R = m_radius;
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
        void get_tail_connections(array_type &col_ptr, array_type &col_idx)
        {
            for(size_t j = 0; j < m_dims[1]; ++j)
                for(size_t i = 0; i < m_dims[0]; ++i)
                {
                    this->add_plane_connections(i, j, m_dims[0], m_dims[1], col_idx);
                    this->add_cylinder_connections(i, j, m_dims[0], m_dims[1], col_idx);
                    col_ptr.push_back(col_idx.size());
                }
        }

        template<typename array_type>
        void get_head_connections(array_type &col_ptr, array_type &col_idx)
        {
            for(size_t i = 0; i < m_dims[2]; ++i)
                col_idx.push_back(i + 1);
            col_ptr.push_back(col_idx.size());
            for(size_t j = 0; j < m_dims[3] - 1; ++j)
                for(size_t i = 0; i < m_dims[2]; ++i)
                {
                    this->add_plane_connections(i, j, m_dims[2], m_dims[3] - 1, col_idx, 1);
                    this->add_cylinder_connections(i, j, m_dims[2], m_dims[3] - 1, col_idx, 1);
                    if(j == m_dims[3] - 2)
                    {
                        size_t head_offset = m_dims[2]*(m_dims[3]-1)+1;
                        for(size_t k = head_offset; k < head_offset + m_dims[0]; ++k)
                            col_idx.push_back(k);
                    }
                    col_ptr.push_back(col_idx.size());
                }            
        }

        void setCells(size_t offset)
        {
            // Set cells for the nose of the head
            {
                for(size_t i = 1; i < m_dims[2]; ++i)
                {
                    vtkIdType cell[3] = {0+offset, i + 1+offset, i+offset};
                    this->getCells()->InsertNextCell(3, cell);
                }
                vtkIdType cell[3] = {0+offset, 1+offset, m_dims[2]+offset};
                this->getCells()->InsertNextCell(3, cell);
            }
            for(size_t j = 0; j < m_dims[3] - 1; ++j)
                for(size_t i = 0; i < m_dims[2]; ++i)
                {
                    this->set_corner_cells(i, j,  m_dims[2], m_dims[3] - 1, 1+offset);
                    this->set_plane_cells(i, j, m_dims[2], m_dims[3] - 1, 1+offset);
                }
            size_t head_offset = m_dims[2] * (m_dims[3] - 1) + 1;
            for(size_t j = 0; j < m_dims[1]; ++j)
                for(size_t i = 0; i < m_dims[0]; ++i)
                {
                    this->set_corner_cells(i, j, m_dims[0], m_dims[1], head_offset+offset);
                    this->set_plane_cells(i, j, m_dims[0], m_dims[1], head_offset+offset);
                }
            setJunctionCells(offset);
        }

    private:

        void setJunctionCells(size_t offset)
        {
            int factor = m_dims[2] / m_dims[0];
            int head_offset = m_dims[2] * (m_dims[3] - 1) + 1;
            int head_points = head_offset - m_dims[2];
            for(size_t i = 0; i < m_dims[0]; ++i)
            {
                vtkIdType cells[3][3] = {{i + head_offset+offset, head_points + i *factor + 1+offset, head_points + i *factor+offset},
                {i + head_offset+offset, (i + 1) % m_dims[0] + head_offset+offset, head_points + i *factor + 1+offset},
                {(i + 1) % m_dims[0] + head_offset+offset, head_points + (i *factor + 2) % m_dims[2]+offset, head_points + i *factor + 1+offset}
                };
                for(int k = 0; k < 3; ++k)
                    this->getCells()->InsertNextCell(3, cells[k]);
            }
        }

        value_type amplitude_cubic(const value_type &s)
        {
            value_type value = -m_amplitude * 0.03125 * (-70 + 84 * s - 35 * std::pow(s, 2) + 5 * std::pow(s, 3)) * std::pow(s, 4);
            return s >= 2 ? m_amplitude : s <= 0 ? 0. : value;
        }

        value_type damplitude_cubic(const value_type &s)
        {
            value_type value = -m_amplitude * 0.03125 * 35 * std::pow(s, 3) * std::pow(s - 2, 3);;
            return (s >= 2 || s <= 0) ? 0. : value;
        }

        void get_tail_frame(value_type s, value_type t, value_type x[2], value_type normal[2])
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


