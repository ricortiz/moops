#ifndef OVAL_GEOMETRY_HPP
#define OVAL_GEOMETRY_HPP

#include <cmath>
#include "base_geometry.hpp"
#include "particle_system/particle.hpp"

template<typename value_type>
class OvalGeometry : public BaseGeometry<OvalGeometry<value_type> >
{

protected:

    typedef BaseGeometry<OvalGeometry<value_type> > base_type;
    typedef Particle<value_type> particle_type;
    
    private:
        value_type       m_x0[3];
        size_t           m_dims[2];
        value_type       m_speed;
        size_t           m_size;
        value_type       m_inner_radius;
        value_type       m_outer_radius;
        size_t           m_hi;
        size_t           m_lo;
        std::vector<value_type> m_radius_scale;

    public:

        OvalGeometry() : m_speed(1), m_size(0), m_inner_radius(.05), m_outer_radius(.5), m_lo(250), m_hi(750) {}
        OvalGeometry(value_type x0[3], size_t M = 6, size_t N = 100, value_type speed = .05, value_type inner_radius = .05, value_type outer_radius = .5)
                : m_speed(speed), m_size(M*N), m_inner_radius(inner_radius), m_outer_radius(outer_radius), m_hi(500), m_lo(1000)
        {
            m_dims[0] = M;
            m_dims[1] = N;
            m_x0[0] = x0[0];
            m_x0[1] = x0[1];
            m_x0[2] = x0[2];
            assert(m_hi > m_lo);
            m_radius_scale.resize(m_hi-m_lo,0.0);
        }

        void getDimensions(size_t &M, size_t &N) { M = m_dims[0]; N = m_dims[1]; }
        void getInnerRadius(value_type &radius) { radius = m_inner_radius; }
        void getOuterRadius(value_type &radius) { radius = m_outer_radius; }
        void getForcingRange(size_t &lo, size_t &hi) { hi = m_hi; lo = m_lo; }
        void getWaveSpeed(value_type &speed) { speed = m_speed; }
        
        void getX0(value_type *x0) { x0 = m_x0; }
        void setDimensions(size_t dims[2]) { m_dims[0] = dims[0]; m_dims[1] = dims[1]; }
        void setDimensions(size_t M, size_t N) { m_dims[0] = M; m_dims[1] = N; }
        void setInnerRadius(value_type radius) { m_inner_radius = radius; }
        void setOuterRadius(value_type radius) { radius = m_outer_radius; }
        void setForcingRange(size_t lo, size_t hi) { m_hi = hi; m_lo = lo; m_radius_scale.resize(m_hi-m_lo,0.0); }
        void setWaveSpeed(value_type speed) { m_speed = speed; }
        void setX0(value_type x[3]) { m_x0[0] = x[0]; m_x0[1] = x[1]; m_x0[2] = x[2]; }
        void setX0(value_type x, value_type y, value_type z)
        {
            m_x0[0] = x;
            m_x0[1] = y;
            m_x0[2] = z;
        }
        
        size_t const &size() const { return m_size; }
        void getLocalFrame(value_type s, value_type x[3],value_type normal[3])
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

        inline void surface_point(size_t i, size_t j, value_type scale, value_type *position, value_type dtheta, value_type dalpha)
        {

            value_type s = i*dalpha;

            /*centerline and tangent*/
            value_type x[3] = {0}, normal[3] = {0};
            getLocalFrame(s,x,normal);

            value_type theta = j * dtheta;
            value_type R = scale*m_inner_radius;
            position[0] = R * normal[0] * std::cos(theta) + x[0] + m_x0[0];
            position[1] = R * normal[1] * std::cos(theta) + x[1] + m_x0[1];
            position[2] = R * std::sin(theta)                    + m_x0[2];
        }
        
        void surface_point(size_t i, size_t j, value_type scale, particle_type &particle, value_type dtheta, value_type dalpha)
        {
           surface_point(i,j,scale,particle.position,dtheta,dalpha);
        }
        
        inline void setCells()
        {
            for (size_t i = 0; i < m_dims[0]; ++i)
                for (size_t j = 0; j < m_dims[1]; ++j)
                {
                    this->set_corner_cells(i,j,m_dims[0],m_dims[1]);
                    this->set_top_cells(i,j,m_dims[0],m_dims[1]);
                    this->set_plane_cells(i,j,m_dims[0],m_dims[1]);
                }
        }

        template<typename array_type>
        inline void getConnections(array_type &col_ptr, array_type &col_idx)
        {
            for (size_t j = 0; j < m_dims[1]; ++j)
                for (size_t i = 0; i < m_dims[0]; ++i)
                {
                    this->add_cylinder_connections(i,j,m_dims[0],m_dims[1],col_idx);
                    this->add_plane_connections(i,j,m_dims[0],m_dims[1],col_idx);
                    this->add_closed_connections(i,j,m_dims[0],m_dims[1],col_idx);
                    col_ptr.push_back(col_idx.size());
                }
        }

        template<typename spring_type>
        inline void resetRestingLength(spring_type &spring, value_type time)
        {            
            int idx = std::max(std::min(spring->A()->i,m_hi-1),m_lo) - m_lo;
            assert( idx >= 0 && idx < m_radius_scale.size() );
            base_type::resetRestingLength(spring,m_radius_scale[idx]);
        }
               
        const std::vector<value_type> &setRadiusScaling(value_type t)
        {
            logger.startTimer("radiusScale");
            size_t n = m_hi - m_lo;
            // Range of stretchy part
            value_type x_1 = 0.0;
            value_type x_2 = 1.0;
            
            // Range of part to be squeezed
            value_type x_a = x_1;
            value_type x_b = x_2;
            
            // size of squeeze
            value_type Ls = 0.3;
            value_type radius = .015;
            
            value_type freq = 3.0;
            value_type wt = std::fmod(t / freq, value_type(1));
            value_type x_c = (1 - wt) * x_a + wt * (x_b - Ls);
            value_type x_d = x_c + Ls;
            
            value_type dx = (x_2 - x_1) / n;
            
            std::vector<value_type> x(n), filter(n);
            value_type scale = 0;
            for(size_t k = 0; k < n; ++k)
            {
                x[k] = dx * k;
                filter[k] = (x[k] > x_c) && (x[k] < x_d);
                m_radius_scale[k] = (x[k] - x_c) * (x[k] - x_d);
                m_radius_scale[k] *= m_radius_scale[k];
                if(filter[k] && (m_radius_scale[k] > scale))
                    scale = m_radius_scale[k];
            }
            scale = std::sin(M_PI * wt) * radius / (scale * m_inner_radius);
            for(size_t k = 0; k < n; ++k)
            {
                m_radius_scale[k] *= filter[k] * scale;
                m_radius_scale[k] = 1 - m_radius_scale[k];
            }
            logger.stopTimer("radiusScale");
            return m_radius_scale;
        }

};

template<typename _value_type>
struct geometry_traits<OvalGeometry<_value_type> >
{
    typedef _value_type          value_type;
    typedef Particle<value_type> particle_type;
};

#endif
