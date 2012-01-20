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

        template<typename particle_type>
        void init_surface(sperm_type &geometry, particle_type *particles, value_type num_particles)
        {
            size_t *dimensions = geometry.get_dimensions();
            size_t points_per_geometry = dimensions[0] * dimensions[1] + dimensions[2] * (dimensions[3] - 1) + 1;
            m_num_geometries = num_particles / points_per_geometry;
            std::vector<value_type> mesh2d;
            setGeometryGrid(mesh2d);
            m_tail_col_ptr.push_back(0);
            m_head_col_ptr.push_back(0);
            grid_type &grid = this->grid();
            for(size_t i = 0, idx = 0; i < m_num_geometries; ++i)
            {
                geometry.init(&particles[idx], grid);
                value_type *translation_point = &mesh2d[3 * i];
                for(size_t k = 0; k < points_per_geometry; ++k, ++idx)
                {
                    particles[idx].position[0] += translation_point[0];
                    particles[idx].position[1] += translation_point[1];
                    particles[idx].position[2] += translation_point[2];
                }
            }
            geometry.get_head_connections(m_head_col_ptr, m_head_col_idx);
            geometry.get_tail_connections(m_tail_col_ptr, m_tail_col_idx);
            m_sperm = &geometry;
	    
        }



//         template<typename geometry_type>
//         void init_volume(BaseGeometry<geometry_type> &geometry, value_type *positions)
//         {
//             geometry.init(positions);
//         }

        template<typename spring_system_type>
        void set_springs(spring_system_type &spring_system)
        {
            typedef typename spring_system_type::particle_type particle_type;
            particle_type *particles = spring_system.particles();
            size_t *dimensions = m_sperm->get_dimensions();
            size_t tail_offset = dimensions[0] * dimensions[1];
            size_t head_offset = dimensions[2] * (dimensions[3] - 1) + 1;
            size_t offset = tail_offset + head_offset;
//             std::cout << "idx = [";
//             for(size_t p = 0; p < m_tail_col_ptr.size() - 1; ++p)                
//                 for(size_t i = m_tail_col_ptr[p], end = m_tail_col_ptr[p + 1]; i < end; ++i)
//                     std::cout << "[" << p+1 << "," << m_tail_col_idx[i]+1 << "];";
//             std::cout << "];" << std::endl;

//             std::cout << "springs = [";
            for(size_t i = 0; i < m_num_geometries; ++i)
            {
                size_t idx = i * offset;
                for(size_t p = 0; p < m_head_col_ptr.size() - 1; ++p)
                    for(size_t i = m_head_col_ptr[p], end = m_head_col_ptr[p + 1]; i < end; ++i)
                        if(!spring_system.exist_spring(&particles[p + idx], &particles[m_head_col_idx[i] + idx]))
                        {
                            spring_system.add_spring(&particles[p + idx], &particles[m_head_col_idx[i] + idx],.1);
//                             std::cout << p + idx + 1 << "," << m_head_col_idx[i] + idx + 1 << ";";
                        }

                for(size_t p = 0; p < m_tail_col_ptr.size() - 1; ++p)
                    for(size_t i = m_tail_col_ptr[p], end = m_tail_col_ptr[p + 1]; i < end; ++i)
                        if(!spring_system.exist_spring(&particles[p + idx + head_offset], &particles[m_tail_col_idx[i] + idx + head_offset]))
                        {
                            spring_system.add_spring(&particles[p + idx + head_offset], &particles[m_tail_col_idx[i] + idx + head_offset],.1,true);
//                             std::cout << p + idx + head_offset+1 << "," << m_tail_col_idx[i] + idx + head_offset+1  << ";";
                        }
            }
            spring_system.fluid_solver().initMaps(spring_system);
//             std::cout << "];\n";
            std::cout << "Created " << spring_system.springs_size() << " springs." << std::endl;
        }

        template<typename spring_system_type>
        void update(spring_system_type &spring_system, value_type time)
        {
            typedef typename spring_system_type::spring_ptr_container spring_ptr_type;
            typedef typename spring_system_type::spring_ptr_container::iterator iterator;
            typedef typename spring_system_type::spring_type spring_type;
            spring_ptr_type &springs_map = spring_system.springs_update_map();
            spring_type *spring;
            for(iterator s = springs_map.begin(), end = springs_map.end(); s != end; ++s)
            {
                spring = *s;
                m_sperm->resetRestingLength(spring,time);
            }
            updateForceGradient(spring_system);
        }


        void setGeometryGrid(std::vector<value_type> &mesh2d)
        {
            value_type dtheta = 1., dalpha = 1.;
            for(size_t i = 0, end = std::sqrt(m_num_geometries); i < end; ++i)
                for(size_t j = 0; j < end; ++j)
                {
                    mesh2d.push_back(i * dtheta);
                    mesh2d.push_back(0.0);
                    mesh2d.push_back(j * dalpha);
                }
        }

        template<typename spring_system_type>
        void updateForceGradient(spring_system_type &spring_system)
        {
            typedef typename spring_system_type::particle_type particle_type;
            size_t *dimensions = m_sperm->get_dimensions();
            size_t tail_offset = dimensions[0] * dimensions[1];
            size_t head_offset = dimensions[2] * (dimensions[3] - 1) + 1;
            size_t offset = tail_offset + head_offset;

            value_type gradient[3] = {0};
            particle_type *particles = spring_system.particles();
            for(size_t i = 0; i < head_offset; ++i)
            {
                gradient[0] += particles[i].force[0];
                gradient[1] += particles[i].force[1];
                gradient[2] += particles[i].force[2];
            }
            gradient[0] /= head_offset;
            gradient[1] /= head_offset;
            gradient[2] /= head_offset;
            value_type norm = std::sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1] + gradient[2] * gradient[2]);
            if(norm == 0) return;
            
            value_type scale = 2.0/* / norm*/;

            gradient[0] = -particles[0].position[0] * scale;
            gradient[1] = -particles[0].position[1] * scale;
            gradient[2] = -particles[0].position[2] * scale;
            particles[0].force[0] += gradient[0];
            particles[0].force[1] += gradient[1];
            particles[0].force[2] += gradient[2];
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
