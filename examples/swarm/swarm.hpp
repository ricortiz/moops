#ifndef SWARM_HPP
#define SWARM_HPP

#include <cassert>
#include <deque>
#include <cmath>
#include "geometry/surface.hpp"
#include "geometry/sine_geometry.hpp"

template<typename value_type>
class Swarm : public Surface<Swarm<value_type> >
{
    public:
        typedef SineGeometry<value_type> sperm_type;
        typedef typename Surface<Swarm<value_type> >::grid_type grid_type;

    private:
        sperm_type *m_sperm;
        std::vector<size_t>    m_tail_col_ptr;
        std::vector<size_t>    m_tail_col_idx;
        std::vector<size_t>    m_head_col_ptr;
        std::vector<size_t>    m_head_col_idx;
        size_t m_num_geometries;

    public:

        void init_surface ( sperm_type &geometry, value_type *positions, value_type num_particles )
        {
            size_t *dimensions = geometry.get_dimensions();
            size_t points_per_geometry = dimensions[0] * dimensions[1] + dimensions[2] * (dimensions[3]-1) +1;
            m_num_geometries = num_particles / points_per_geometry;
            m_tail_col_ptr.push_back ( 0 );
            m_head_col_ptr.push_back ( 0 );
            grid_type &grid = this->grid();
            for ( size_t i = 0; i < m_num_geometries; ++i )
                geometry.init ( positions + i*points_per_geometry, grid );
//             std::cout << "p = [";
//             for(size_t i = 0; i < 3*points_per_geometry-1; ++i)
//                 std::cout << positions[i] << ",";
//             std::cout << positions[3*points_per_geometry-1] << "];" << std::endl;
            geometry.get_head_connections ( m_head_col_ptr, m_head_col_idx );
            geometry.get_tail_connections ( m_tail_col_ptr, m_tail_col_idx );
            m_sperm = &geometry;
        }

        template<typename geometry_type>
        void init_volume ( BaseGeometry<geometry_type> &geometry, value_type *positions )
        {
            geometry.init ( positions );
        }

        template<typename spring_system_type>
        void set_springs ( spring_system_type &spring_system )
        {
            typedef typename spring_system_type::particle_type particle_type;

            particle_type *particles = spring_system.particles();
            size_t *dimensions = m_sperm->get_dimensions();
            size_t tail_offset = dimensions[0] * dimensions[1];
            size_t head_offset = dimensions[2] * (dimensions[3]-1) +1;
            size_t offset = tail_offset + head_offset;

            std::cout << "springs = [";
            for ( size_t i = 0; i < m_num_geometries; ++i )
            {
                size_t idx = i * offset;
                for ( size_t p = 0; p < m_head_col_ptr.size() - 1; ++p )
                    for ( size_t i = m_head_col_ptr[p], end = m_head_col_ptr[p+1]; i < end; ++i )
                        if ( !spring_system.exist_spring ( &particles[p+idx], &particles[m_head_col_idx[i]+idx] ) )
                            if (spring_system.add_spring ( &particles[p+idx], &particles[m_head_col_idx[i]+idx], 1. ))
                                std::cout << p+1 << "," << m_head_col_idx[i]+1 << ";";
                std::cout << "];\n";
//                 for ( size_t p = 0; p < m_tail_col_ptr.size() - 1; ++p )
//                     for ( size_t i = m_tail_col_ptr[p], end = m_tail_col_ptr[p+1]; i < end; ++i )
//                         if ( !spring_system.exist_spring ( &particles[p+idx+head_offset], &particles[m_tail_col_idx[i]+idx+head_offset] ) )
//                             spring_system.add_spring ( &particles[p+idx+head_offset], &particles[m_tail_col_idx[i]+idx+head_offset], 1. );
            }
            std::cout << "Created " << spring_system.springs_size() << " springs." << std::endl;
        }

        template<typename spring_system_type>
        void update ( spring_system_type &spring_system, value_type time )
        {
            typedef typename spring_system_type::spring_iterator iterator;
            grid_type &grid = this->grid();

            for ( iterator s = spring_system.springs_begin(), end = spring_system.springs_end(); s != end; ++s )
            {
                size_t Ai = grid[ s->A()->position].first;
                size_t Aj = grid[ s->A()->position].second;
                size_t Bi = grid[ s->B()->position].first;
                size_t Bj = grid[ s->B()->position].second;
                s->resting_length() = m_sperm->get_distance ( Ai, Aj, Bi, Bj, time );
            }

        }

        inline BaseGeometry<sperm_type> *geometry()
        {
            return m_sperm;
        }

};


template<typename _value_type>
struct surface_traits<Swarm<_value_type> >
{
    typedef _value_type value_type;
    typedef SineGeometry<value_type> geometry_type;
};

#endif
