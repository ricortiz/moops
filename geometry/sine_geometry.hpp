#ifndef SINE_GEOMETRY_HPP
#define SINE_GEOMETRY_HPP

#include<vector>
#include <cmath>
#include "base_geometry.hpp"
#include <boost/type_traits/detail/is_mem_fun_pointer_impl.hpp>

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

        SineGeometry ( )
            : m_length ( 4.0 ), m_speed ( .5 ), m_amplitude ( .25 ), m_pitch ( 4.1 ), m_radius ( .05 ), m_size ( 0 ) {}
        SineGeometry ( value_type x0[3], size_t M = 6, size_t N = 100, value_type length = 4.0, value_type speed = .5, value_type amplitude = .25, value_type pitch = 4.1, value_type radius = .05 )
            : m_length ( length ), m_speed ( speed ), m_amplitude ( amplitude ), m_pitch ( pitch ), m_radius ( radius ), m_size ( 0 )
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
        value_type *get_x0() { return m_x0; }

        void set_x0 ( value_type x, value_type y, value_type z )
        {
            m_x0[0] = x;
            m_x0[1] = y;
            m_x0[2] = z;
        }

        void surface_point ( size_t i, size_t j, value_type t, value_type *point, value_type dtheta, value_type dalpha )
        {
            value_type s = i * dalpha;
            /*centerline and tangent*/
            value_type x[2] = {0}, normal[2] = {0};
            get_tail_frame ( s, t, x, normal );
            value_type theta = j * dtheta;
            point[0] = m_radius * normal[0] * std::cos ( theta ) + x[0] + m_x0[0];
            point[1] = m_radius * normal[1] * std::cos ( theta ) + x[1] + m_x0[1];
            point[2] = m_radius * std::sin ( theta )                    + m_x0[2];
            m_size++;
        }

        template<typename array_type>
        void get_tail_connections ( array_type &col_ptr, array_type &col_idx )
        {
            for ( size_t j = 0; j < m_dims[1]; ++j )
                for ( size_t i = 0; i < m_dims[0]; ++i )
                {
                    this->add_plane_connections ( i, j, m_dims[0], m_dims[1], col_idx );
                    this->add_cylinder_connections ( i, j, m_dims[0], m_dims[1], col_idx );
                    col_ptr.push_back ( col_idx.size() );
                }
        }

        template<typename array_type>
        void get_head_connections ( array_type &col_ptr, array_type &col_idx )
        {
            for ( size_t i = 0; i < m_dims[2]; ++i )
                col_idx.push_back(i+1);
            col_ptr.push_back ( col_idx.size());
            for ( size_t j = 0; j < 1; ++j )
            {
                for ( size_t i = 0; i < m_dims[2]; ++i )
                {
                    this->add_plane_connections ( i, j, m_dims[2], m_dims[3]-1, col_idx, 1 );
                    this->add_cylinder_connections ( i, j, m_dims[2], m_dims[3]-1, col_idx, 1 );
                    col_ptr.push_back ( col_idx.size() );
                }
            }
        }

        void head_point ( size_t i, size_t j, value_type t, value_type *point, value_type dtheta, value_type dalpha )
        {
            value_type scale = 5.0 * m_radius;
            value_type R = i < 7 ? scale * ( std::cos ( ( 7 - i ) * .5 * M_PI / 7 ) ) : i > 13 ? scale * ( std::cos ( ( i - 13 ) * .5 * M_PI / 6 ) ) : scale;
            if ( i == m_dims[3] - 2 || i == m_dims[3] - 1 )
                R = m_radius;
            value_type s = ( i - .8 ) * dalpha;
            /*centerline and tangent*/

            value_type normal = 1.0;
            value_type binormal = 1.0;
            value_type theta = j * dtheta;
            point[0] = R * normal * std::cos ( theta )     + m_x0[0];
            point[1] =                                   s + m_x0[1];
            point[2] = R * binormal * std::sin ( theta )   + m_x0[2];
            m_size++;
        }

        template<typename grid_type>
        void init ( value_type *positions, grid_type &grid )
        {
            value_type dtheta = 2 * M_PI / m_dims[2];
            value_type dalpha = m_length / 4.0 / m_dims[3];
            size_t idx = 0;
            positions[idx++] = 0;
            positions[idx++] = 0;
            positions[idx++] = 0;
            grid[&positions[0]] = std::make_pair ( 0, 0 );
            for ( size_t i = 1; i < m_dims[3]; ++i )
                for ( size_t j = 0; j < m_dims[2]; ++j, idx += 3 )
                {
                    head_point ( i, j, 0.0, &positions[idx], dtheta, dalpha );
                    grid[&positions[idx]] = std::make_pair ( i, j );
                }
            m_x0[1] = m_length / 4.0 - .3 * dalpha;
            dtheta = 2 * M_PI / m_dims[0];
            dalpha = 2 * M_PI / m_dims[1];
            for ( size_t i = 0; i < m_dims[1]; ++i )
                for ( size_t j = 0; j < m_dims[0]; ++j, idx += 3 )
                {
                    surface_point ( i, j, 0.0, &positions[idx], dtheta, dalpha );
                    grid[&positions[idx]] = std::make_pair ( i, j );
                }
        }

        void setCells()
        {
            for ( size_t i = 0; i < m_dims[0]; ++i )
                for ( size_t j = 0; j < m_dims[1]; ++j )
                {
                    this->set_corner_cells ( i, j, m_dims[0], m_dims[1] );
                    this->set_plane_cells ( i, j, m_dims[0], m_dims[1] );
                }
            size_t offset = m_dims[0] * m_dims[1];
            for ( size_t i = 0; i < m_dims[2]; ++i )
                for ( size_t j = 0; j < m_dims[3]; ++j )
                {
                    this->set_corner_cells ( i, j, m_dims[0], m_dims[1], offset );
                    this->set_plane_cells ( i, j, m_dims[0], m_dims[1], offset );
                }
        }

    private:

        value_type amplitude_cubic ( const value_type &s )
        {
            value_type value = -m_amplitude * 0.03125 * ( -70 + 84 * s - 35 * std::pow ( s, 2 ) + 5 * std::pow ( s, 3 ) ) * std::pow ( s, 4 );
            return s >= 2 ? m_amplitude : s <= 0 ? 0. : value;
        }

        value_type damplitude_cubic ( const value_type &s )
        {
            value_type value = -m_amplitude * 0.03125 * 35 * std::pow ( s, 3 ) * std::pow ( s - 2, 3 );;
            return ( s >= 2 || s <= 0 ) ? 0. : value;
        }

        void get_tail_frame ( value_type s, value_type t, value_type x[2], value_type normal[2] )
        {
            value_type A = amplitude_cubic ( s );
            value_type DA = damplitude_cubic ( s );
            value_type Dx[2];

            x[0] = A * std::cos ( m_pitch * s + t * m_speed );
            x[1] = s;

            /*tangential vector: derivatives*/
            Dx[0] = DA * std::cos ( m_pitch * s + t * m_speed ) - m_pitch * A * std::sin ( m_pitch * s + t * m_speed );
            Dx[1] = 1;

            value_type norm = std::sqrt ( Dx[0] * Dx[0] + Dx[1] * Dx[1] );
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


