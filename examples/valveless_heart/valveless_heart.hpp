#ifndef HEART_PUMP_HPP
#define HEART_PUMP_HPP

#include <cassert>
#include <deque>
#include <cmath>
#include <iterator>
#include "geometry/surface.hpp"
#include "geometry/torus_geometry.hpp"
#include "geometry/oval_geometry.hpp"


template<typename value_type>
class HeartPump : public Surface<HeartPump<value_type> >
{
public:
    typedef TorusGeometry<value_type> torus_type;
    typedef OvalGeometry<value_type> oval_type;
    typedef typename Surface<HeartPump<value_type> >::grid_type grid_type;

private:
    oval_type *m_geometry;
    std::vector<size_t>    m_col_ptr;
    std::vector<size_t>    m_col_idx;
    size_t m_hi;
    size_t m_lo;
    std::vector<value_type> m_radius_scale;

public:

    template<typename particle_type>
    void init_surface(oval_type &oval_geometry, particle_type *particles)
    {
        m_lo = 10; m_hi = 60;
        m_radius_scale.resize(m_hi-m_lo);
        m_col_ptr.push_back(0);
        grid_type &grid = this->grid();
        oval_geometry.init(particles,grid);
        oval_geometry.get_connections(m_col_ptr,m_col_idx);
        m_geometry = &oval_geometry;
    }
    template<typename particle_type>
    void init_volume(BaseGeometry<oval_type> &oval_geometry, particle_type *particles, size_t num_sub_surfaces = 1)
    {
        oval_geometry.init(particles,num_sub_surfaces);
    }

    template<typename spring_system_type>
    void set_springs(spring_system_type &spring_system)
    {
        typedef typename spring_system_type::particle_type particle_type;
        particle_type *particles = spring_system.particles();

        // Use this to determine wich cross-sections of the geometry are going to be stiffer
        size_t lo = (m_lo+2)*20, hi = (m_hi-2)*20;
//         std::cout << "springs = [";
        for (size_t p = 0; p < m_col_ptr.size()-1; ++p)
        {
            for (size_t i = m_col_ptr[p], end = m_col_ptr[p+1]; i < end; ++i)
                if (!spring_system.exist_spring(&particles[p],&particles[m_col_idx[i]]))
                {
                    if (p > lo && p < hi){
                        /*if (*/spring_system.add_spring(&particles[p],&particles[m_col_idx[i]],2.)/*)
                            std::cout << p+1 << "," << m_col_idx[i]+1 << ";"*/;
                    }
                    else{
                        /*if (*/spring_system.add_spring(&particles[p],&particles[m_col_idx[i]],3.)/*)
                            std::cout << p+1 << "," << m_col_idx[i]+1 << ";"*/;
                    }
                    
                }
        }
//         std::cout << "];\n";
        std::cout << "Created " << spring_system.springs_size() << " springs." << std::endl;
    }

    template<typename spring_system_type>
    void update(spring_system_type &spring_system, value_type time)
    {
        typedef typename spring_system_type::spring_ptr_container::iterator iterator;
        typedef typename spring_system_type::spring_lut_type spring_map_type;
        typedef typename spring_system_type::particle_type particle_type;
        grid_type &grid = this->grid();
        spring_map_type &springs_map = spring_system.springs_map();
        particle_type *particles = spring_system.particles();
        size_t *dims = m_geometry->get_dimensions();

        radius_scale(time);
        for(size_t i = m_lo, idx = m_lo*dims[0]; i < m_hi; ++i)
            for(size_t j = 0; j < dims[0]; ++j, ++idx)
            {
                particle_type *particle = &particles[idx];
                for (iterator s = springs_map[particle].begin(), end = springs_map[particle].end(); s != end; ++s)
                {
                    size_t Ai = grid[(*s)->A()].first;
                    size_t Aj = grid[(*s)->A()].second;
                    size_t Bi = grid[(*s)->B()].first;
                    size_t Bj = grid[(*s)->B()].second;
                    (*s)->resting_length() = m_geometry->get_distance(Ai,Aj,Bi,Bj,m_radius_scale[i-m_lo]);
                }
            }
    }

    std::vector<value_type> &radius_scale(value_type t)
    {
        size_t n = m_hi-m_lo;
        // Range of stretchy part
        value_type x_1 = 0.0;
        value_type x_2 = 1.0;
        
        // Range of part to be squeezed
        value_type x_a = x_1;
        value_type x_b = x_2;
        
        // size of squeeze
        value_type Ls = 0.3;
        value_type radius = .035;
        
        value_type freq = 3.0;
        value_type wt = std::fmod(t/freq,value_type(1));
        value_type x_c = (1-wt)*x_a + wt*(x_b-Ls);
        value_type x_d = x_c + Ls;
        
        value_type dx = (x_2-x_1)/n;
        
        std::vector<value_type> x(n), filter(n);
        value_type scale = 0;
        for(size_t k = 0; k < n; ++k)
        {
            x[k] = dx*k;
            filter[k] = (x[k] > x_c) && (x[k] < x_d);
            m_radius_scale[k] = (x[k]-x_c)*(x[k]-x_d);
            m_radius_scale[k] *= m_radius_scale[k];
            if(filter[k] && (m_radius_scale[k] > scale))
                scale = m_radius_scale[k];
        }
        scale = std::sin(M_PI*wt)*radius/(scale*m_geometry->inner_radius());
        for(size_t k = 0; k < n; ++k)
        {            
            m_radius_scale[k] *= filter[k]*scale;
            m_radius_scale[k] = 1 - m_radius_scale[k];
        }
        return m_radius_scale;
    }
    
    inline BaseGeometry<oval_type> *geometry() {
        return m_geometry;
    }

    void get_lohi(size_t &lo,size_t &hi) { lo = m_lo; hi = m_hi; }

};


template<typename _value_type>
struct surface_traits<HeartPump<_value_type> >
{
    typedef _value_type value_type;
    typedef OvalGeometry<value_type> geometry_type;
};

#endif
