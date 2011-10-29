#ifndef TORUS_GEOMETRY_HPP
#define TORUS_GEOMETRY_HPP

#include<vector>
#include <cmath>
#include "base_geometry.hpp"

template<typename value_type>
class TorusGeometry : public BaseGeometry<TorusGeometry<value_type> >
{

    private:
        value_type       m_x0[3];
        size_t           m_dims[2];
        value_type       m_speed;
        size_t           m_size;
        value_type       m_inner_radius;
        value_type       m_outer_radius;

    public:

        TorusGeometry() : m_speed(1), m_inner_radius(.05), m_outer_radius(.5), m_size(0) {}
        TorusGeometry(value_type x0[3], size_t M = 6, size_t N = 100, value_type speed = .0001, value_type inner_radius = .05, value_type outer_radius = .5)
            : m_speed(speed), m_size(M*N), m_inner_radius(inner_radius), m_outer_radius(outer_radius)
        {
            m_dims[0] = M;
            m_dims[1] = N;
            m_x0[0] = x0[0];
            m_x0[1] = x0[1];
            m_x0[2] = x0[2];
        }

        
        size_t *get_dimensions() { return m_dims; }
        size_t const *get_dimensions() const { return m_dims; }
        size_t size() { return m_size; }
        value_type *get_x0() { return m_x0; }

        void set_x0(value_type x, value_type y, value_type z)
        {
            m_x0[0] = x;
            m_x0[1] = y;
            m_x0[2] = z;
        }

        inline void surface_point(size_t i, size_t j, value_type t, value_type *point)
        {
            value_type dtheta = 2 * M_PI / m_dims[0];
            value_type dalpha = 2 * M_PI / m_dims[1];

            value_type s = i*dalpha;

            /*centerline and tangent*/
            value_type x[3] = {0}, Dx[3] = {0};

            x[0] = m_outer_radius * std::cos(s);
            x[1] = m_outer_radius * std::sin(s);

            /*tangential vector: partial derivatives*/
            Dx[0] = -m_outer_radius * std::sin(s);
            Dx[1] = m_outer_radius * std::cos(s);

            value_type norm = std::sqrt(Dx[0]*Dx[0]+Dx[1]*Dx[1]);
            Dx[0]/=norm;
            Dx[1]/=norm;

            /*normal vector*/
            value_type Normal[3] = {Dx[1],-Dx[0]};

            /*binormal vector: cross product of n and t */
            value_type Binormal[3] = {0.0,0.0,Normal[0]*Dx[1]-Normal[1]*Dx[0]};

            value_type theta = j * dtheta;
            value_type R = .5*m_inner_radius * std::cos(m_speed*t) * std::cos(m_speed*t)+.5*m_inner_radius;
            point[0] = R * Normal[0] * std::cos(theta) + x[0] + m_x0[0];
            point[1] = R * Normal[1] * std::cos(theta) + x[1] + m_x0[1];
            point[2] = R * Binormal[2]* std::sin(theta) + m_x0[2];
        }

        inline void volume_point(size_t i, size_t j, value_type *point, value_type depth = .1)
        {
            value_type dtheta = M_PI / m_dims[0];
            value_type dalpha = M_PI / m_dims[1];
            
            value_type s = i*dalpha;
            
            /*centerline and tangent*/
            value_type x[3] = {0}, Dx[3] = {0};
            
            x[0] = m_outer_radius * std::cos(s);
            x[1] = m_outer_radius * std::sin(s);
            
            /*tangential vector: partial derivatives*/
            Dx[0] = -m_outer_radius * std::sin(s);
            Dx[1] = m_outer_radius * std::cos(s);
            
            value_type norm = std::sqrt(Dx[0]*Dx[0]+Dx[1]*Dx[1]);
            Dx[0]/=norm;
            Dx[1]/=norm;
            
            /*normal vector*/
            value_type Normal[3] = {Dx[1],-Dx[0]};
            
            /*binormal vector: cross product of n and t */
            value_type Binormal[3] = {0.0,0.0,Normal[0]*Dx[1]-Normal[1]*Dx[0]};
            
            value_type theta = j * dtheta;
            value_type R = depth*m_inner_radius;
            point[0] = R * Normal[0] * std::cos(theta) + x[0] + m_x0[0];
            point[1] = R * Normal[1] * std::cos(theta) + x[1] + m_x0[1];
            point[2] = R * Binormal[2]* std::sin(theta) + m_x0[2];
        }
        
        template<typename vtk_cell_type>
        inline void get_vtk_cells(vtk_cell_type &cells)
        {
            for(size_t i = 0; i < m_dims[0]; ++i)
                for(size_t j = 0; j < m_dims[1]; ++j)
                {
                    this->set_corner_cells(i,j,cells);
                    this->set_top_cells(i,j,cells);
                    this->set_plane_cells(i,j,cells);
                }
        }

        template<typename array_type>
        inline void get_connections(array_type &col_ptr, array_type &col_idx)
        {
            for (size_t j = 0; j < m_dims[1]; ++j)
                for (size_t i = 0; i < m_dims[0]; ++i)
                {
                    this->add_cylinder_connections(i,j,col_idx);
                    this->add_plane_connections(i,j,col_idx);
                    this->add_closed_connections(i,j,col_idx);
                    col_ptr.push_back(col_idx.size());
                }
        }

};

template<typename _value_type>
struct geometry_traits<TorusGeometry<_value_type> >
{
    typedef _value_type          value_type;
};

#endif
