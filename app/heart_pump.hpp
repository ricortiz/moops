#ifndef HEART_PUMP_HPP
#define HEART_PUMP_HPP

#include <cassert>
#include <deque>
#include <cmath>
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

public:

    void init_surface(oval_type &oval_geometry, value_type *positions)
    {
        m_col_ptr.push_back(0);
        grid_type &grid = this->grid();
        oval_geometry.init(positions,grid);
        oval_geometry.get_connections(m_col_ptr,m_col_idx);
        m_geometry = &oval_geometry;
    }

    void init_volume(BaseGeometry<oval_type> &oval_geometry, value_type *positions, size_t num_sub_surfaces = 1)
    {
        oval_geometry.init(positions,num_sub_surfaces);
    }

    template<typename spring_system_type>
    void set_springs(spring_system_type &spring_system)
    {
        typedef typename spring_system_type::particle_type particle_type;
        particle_type *particles = spring_system.particles();

        // Use this to determine wich cross-sections of the geometry are going to be stiffer
        size_t lo = 400, hi = 1000;
        for (size_t p = 0; p < m_col_ptr.size()-1; ++p)
        {
            for (size_t i = m_col_ptr[p], end = m_col_ptr[p+1]; i < end; ++i)
                if (!spring_system.exist_spring(&particles[p],&particles[m_col_idx[i]]))
                    if (p > lo && p < hi)
                        spring_system.add_spring(&particles[p],&particles[m_col_idx[i]],1.);
                    else
                        spring_system.add_spring(&particles[p],&particles[m_col_idx[i]],10.);
        }
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
        const size_t lo = 20, hi = 50;
        size_t *dims = m_geometry->get_dimensions();

        std::vector<value_type> rscale(hi-lo, 0.0);
        radius_scale(time,lo,hi,rscale);
        for(size_t i = lo, idx = lo*dims[0]; i < hi; ++i)
            for(size_t j = 0; j < dims[0]; ++j, ++idx)
            {
                particle_type *particle = &particles[idx];
                for (iterator s = springs_map[particle].begin(), end = springs_map[particle].end(); s != end; ++s)
                {
                    size_t Ai = grid[(*s)->A()->position].first;
                    size_t Aj = grid[(*s)->A()->position].second;
                    size_t Bi = grid[(*s)->B()->position].first;
                    size_t Bj = grid[(*s)->B()->position].second;
                    (*s)->resting_length() = m_geometry->get_distance(Ai,Aj,Bi,Bj,rscale[i-lo]);
                }
            }
    }

    void radius_scale(value_type t, size_t lo, size_t hi, std::vector<value_type> &r)
    {        
        // Range of stretchy part
        value_type x_1 = 0.0;
        value_type x_2 = 1.0;
        
        // Range of part to be squeezed
        value_type x_a = 0.2;
        value_type x_b = 0.8;
        
        // size of squeeze
        value_type Ls = 0.3;
        value_type radius = .035;
        value_type deltar = m_geometry->inner_radius() - radius;
        
        value_type freq = 2.0;
        value_type wt = std::fmod(t/freq,value_type(1));
        value_type x_c = (1-wt)*x_a + wt*(x_b-Ls);
        value_type x_d = x_c + Ls;
        
        value_type dx = (x_2-x_1)/(hi-lo);
        
        std::vector<value_type> x(hi-lo), q(hi-lo), filter(hi-lo);
        for(size_t k = 0; k < hi-lo; ++k)
        {
            x[k] = dx*k;
            filter[k] = (x[k] > x_c) && (x[k] < x_d);
            q[k] = (x[k]-x_c)*(x[k]-x_d);
            q[k] *= q[k];
        }
        value_type scale = 0;
        for(size_t k = 0; k < hi-lo; ++k)
        {
            if(filter[k])
                if(q[k] > scale)
                    scale = q[k];
        }
        for(size_t k = 0; k < hi-lo; ++k)
        {
            q[k] /= scale;
            r[k] = 1 - std::sin(M_PI*wt)*deltar*q[k]*filter[k]/m_geometry->inner_radius();
        }
        
    }
    
    inline BaseGeometry<oval_type> *geometry() {
        return m_geometry;
    }

};


template<typename _value_type>
struct surface_traits<HeartPump<_value_type> >
{
    typedef _value_type value_type;
    typedef OvalGeometry<value_type> geometry_type;
};

#endif
