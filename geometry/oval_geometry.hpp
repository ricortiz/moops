#ifndef OVAL_GEOMETRY_HPP
#define OVAL_GEOMETRY_HPP

#include <cmath>
#include "base_geometry.hpp"

template<typename value_type>
class OvalGeometry : public BaseGeometry<OvalGeometry<value_type> >
{

    private:
        value_type       m_x0[3];
        size_t           m_dims[2];
        value_type       m_speed;
        size_t           m_size;
        value_type       m_inner_radius;
        value_type       m_outer_radius;
        size_t           m_hi;
        size_t           m_lo;

    public:

        OvalGeometry() : m_speed(1), m_size(0), m_inner_radius(.05), m_outer_radius(.5), m_hi(250), m_lo(750) {}
        OvalGeometry(value_type x0[3], size_t M = 6, size_t N = 100, value_type speed = .0001, value_type inner_radius = .05, value_type outer_radius = .5)
                : m_speed(speed), m_size(M*N), m_inner_radius(inner_radius), m_outer_radius(outer_radius), m_hi(500), m_lo(1000)
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
        value_type &inner_radius() { return m_inner_radius; }
        size_t const &inner_radius() const { return m_inner_radius; }
        value_type &outer_radius() { return m_outer_radius; }
        size_t const &outer_radius() const { return m_outer_radius; }

        void set_x0(value_type x, value_type y, value_type z)
        {
            m_x0[0] = x;
            m_x0[1] = y;
            m_x0[2] = z;
        }

        inline void get_local_frame(value_type s, value_type x[3],value_type normal[3])
        {
            value_type Dx[3] = {0};
            value_type b = m_outer_radius;
            value_type a = b/2;
            value_type a2 = a*a;
            value_type a4 = a2*a2;
            value_type b4 = b*b*b*b;
            value_type c =a2*std::cos(2*s);
            value_type d = std::sqrt((-a4+b4)+a4*(std::cos(2*s))*(std::cos(2*s)));
            value_type M = c+d;
            value_type DM = -2*a2*std::sin(2*s)*(1+c/d);
            value_type sqrtM = std::sqrt(M);

            x[0] = sqrtM * std::cos(s);
            x[1] = sqrtM * std::sin(s);

            /*tangential vector: derivatives*/
            value_type DsqrtM = .5*DM/sqrtM;
            Dx[0] = DsqrtM*std::cos(s)-sqrtM*std::sin(s);
            Dx[1] = DsqrtM*std::sin(s)+sqrtM*std::cos(s);

            value_type norm = std::sqrt(Dx[0]*Dx[0]+Dx[1]*Dx[1]);
            Dx[0]/=norm;
            Dx[1]/=norm;

            /*normal vector*/
            normal[0] = Dx[1];
            normal[1] =-Dx[0];
        }

        inline void surface_point(size_t i, size_t j, value_type scale, value_type *point, value_type dtheta, value_type dalpha)
        {

            value_type s = i*dalpha;

            /*centerline and tangent*/
            value_type x[3] = {0}, normal[3] = {0};
            get_local_frame(s,x,normal);

            value_type theta = j * dtheta;
            value_type R = scale*m_inner_radius;
            point[0] = R * normal[0] * std::cos(theta) + x[0] + m_x0[0];
            point[1] = R * normal[1] * std::cos(theta) + x[1] + m_x0[1];
            point[2] = R * std::sin(theta)                    + m_x0[2];
        }

        inline void setCells()
        {
            for (size_t i = 0; i < m_dims[0]; ++i)
                for (size_t j = 0; j < m_dims[1]; ++j)
                {
                    this->set_corner_cells(i,j);
                    this->set_top_cells(i,j);
                    this->set_plane_cells(i,j);
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

        inline void get_forcing_range(size_t &hi, size_t &lo) { hi = m_hi; lo = m_lo; }

};

template<typename _value_type>
struct geometry_traits<OvalGeometry<_value_type> >
{
    typedef _value_type          value_type;
};

#endif
