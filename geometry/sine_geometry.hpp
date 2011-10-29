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
        size_t           m_dims[2];
        size_t           m_size;
        value_type       m_x0[3];

    public:

        
        SineGeometry(value_type x0[3], size_t M = 6, size_t N = 100, value_type length = 4.0, value_type speed = .5, value_type amplitude = .25, value_type pitch = 4.1, value_type radius = .05)
        : m_length(length), m_speed(speed), m_amplitude(amplitude), m_pitch(pitch), m_radius(radius)
        {
            m_dims[0] = M;
            m_dims[1] = N;
            m_size = m_dims[0]*m_dims[1];
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

        void set_x0(value_type x, value_type y, value_type z)
        {
            m_x0[0] = x;
            m_x0[1] = y;
            m_x0[2] = z;
        }
        value_type amplitude_cubic(const value_type &s)
        {
            value_type value = -m_amplitude*0.03125* (-70+84*s-35*std::pow(s,2) +5*std::pow(s,3)) *std::pow(s,4);
            return s >= 2 ? m_amplitude : s <= 0 ? 0. : value;
        }

        value_type damplitude_cubic(const value_type &s)
        {
            value_type value = -m_amplitude*0.03125*35*std::pow(s,3) *std::pow(s-2,3);;
            return (s >= 2 || s <= 0) ? 0. : value;
        }

        void surface_point(size_t i, size_t j, value_type t, value_type *point)
        {
            value_type dtheta = 2 * M_PI / m_dims[0];
            value_type ds = m_length / (m_dims[1] - 1);

            value_type s = i*ds;

            value_type A = amplitude_cubic(s);
            value_type DA = damplitude_cubic(s);

            /*centerline and tangent*/
            value_type x[3] = {0}, Dx[3] = {0};

            x[0] = A * std::cos(m_pitch * s + t*m_speed);
            x[2] = s;

            /*tangential vector: partial derivatives*/
            Dx[0] = DA * std::cos(m_pitch * s + t*m_speed) - m_pitch * A * std::sin(m_pitch * s + t*m_speed);
            Dx[2] = 1.;

            value_type norm = std::sqrt(Dx[0]*Dx[0]+Dx[2]*Dx[2]);
            Dx[0]/=norm;
            Dx[2]/=norm;

            /*normal vector*/
            value_type Normal[3] = {Dx[2],0.0,-Dx[0]};

            /*binormal vector: cross product of n and t */
            value_type Binormal[3] = {0.0,-Normal[0]*Dx[2]+Normal[2]*Dx[0],0.0};

            value_type theta = j * dtheta;
            point[0] = m_radius * Normal[0] * std::cos(theta) + x[0] + m_x0[0];
            point[1] = m_radius * Binormal[1] * std::sin(theta) + m_x0[1];
            point[2] = m_radius * Normal[2] * std::cos(theta) + x[2] + m_x0[2];
        }

};

template<typename _value_type>
struct geometry_traits<SineGeometry<_value_type> >
{
    typedef _value_type          value_type;
};

#endif
