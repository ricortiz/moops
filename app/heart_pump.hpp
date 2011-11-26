#ifndef HEART_PUMP_HPP
#define HEART_PUMP_HPP

#include <cassert>
#include <deque>
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
        BaseGeometry<oval_type> *m_geometry;
        std::vector<size_t>    m_col_ptr;
        std::vector<size_t>    m_col_idx;

    public:

        void init_surface(BaseGeometry<oval_type> &oval_geometry, value_type *positions)
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
            
            size_t lo, hi;
            m_geometry->get_forcing_range(lo,hi);
//             std::cout << "springs = [";
            for (size_t p = 0; p < m_col_ptr.size()-1; ++p)
            {
                for (size_t i = m_col_ptr[p], end = m_col_ptr[p+1]; i < end; ++i)
                    if (!spring_system.exist_spring(&particles[p],&particles[m_col_idx[i]]))
                    {
                        if (p < hi+9 && p > lo-1)
                        {
//                             if(
                            spring_system.add_spring(&particles[p],&particles[m_col_idx[i]],1.)/*)*/;
//                             std::cout << "[" << p+1 << "," << m_col_idx[i]+1 << "];";
                        }
                        else
                        {
                            /*if(*/
                            spring_system.add_spring(&particles[p],&particles[m_col_idx[i]],10.)/*)*/;
//                             std::cout << "[" << p+1 << "," << m_col_idx[i]+1 << "];";
                        }
                    }
            }



            std::cout << "Created " << spring_system.springs_size() << " springs." << std::endl;
        }

        template<typename spring_system_type>
        void update(spring_system_type &spring_system, value_type time)
        {
            typedef typename spring_system_type::spring_ptr_container::iterator iterator;
//             typedef typename spring_system_type::spring_container::iterator iterator;
            typedef typename spring_system_type::spring_lut_type spring_map_type;
            typedef typename spring_system_type::particle_type particle_type;
            grid_type &grid = this->grid();
            spring_map_type &springs_map = spring_system.springs_map();
            particle_type *particles = spring_system.particles();
//             std::cout.precision(16);
//             std::cout << "resting lengths before = ";
//             for (iterator s = spring_system.springs_begin(), end = spring_system.springs_end(); s != end; ++s)
//             {
//                 std::cout << s->resting_length() << ",";
//             }
//             std::cout << std::endl;
            size_t lo, hi;
            m_geometry->get_forcing_range(lo,hi);
            for(size_t i = lo; i <= lo+80; ++i)
            {
                for(iterator s = springs_map[&particles[i]].begin(), end = springs_map[&particles[i]].end(); s != end; ++s)
                {
                    size_t Ai = grid[(*s)->A()->position].first;
                    size_t Aj = grid[(*s)->A()->position].second;
                    size_t Bi = grid[(*s)->B()->position].first;
                    size_t Bj = grid[(*s)->B()->position].second;
                    (*s)->resting_length() = m_geometry->get_distance(Ai,Aj,Bi,Bj,time);
                }
            }
    




            
//             size_t idx = 0;
//             for (iterator s = spring_system.springs_begin(), end = spring_system.springs_end(); s != end; ++s, ++idx)
//             {
//                 if (idx < 3000 && idx > 1000)
//                 {
//                     size_t Ai = grid[s->A()->position].first;
//                     size_t Aj = grid[s->A()->position].second;
//                     size_t Bi = grid[s->B()->position].first;
//                     size_t Bj = grid[s->B()->position].second;
//                     s->resting_length() = m_geometry->get_distance(Ai,Aj,Bi,Bj,time);
//                 }
//             }
//             std::cout << "resting lengths after = ";
//             for (iterator s = spring_system.springs_begin(), end = spring_system.springs_end(); s != end; ++s)
//             {
//                 std::cout << s->resting_length() << ",";
//             }
//             std::cout << std::endl;

        }

        inline BaseGeometry<oval_type> *geometry() { return m_geometry; }

};


template<typename _value_type>
struct surface_traits<HeartPump<_value_type> >
{
    typedef _value_type value_type;
    typedef OvalGeometry<value_type> geometry_type;
};

#endif
