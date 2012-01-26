#ifndef SWARM_HPP
#define SWARM_HPP

#include "particle_system/storage/particle_system_storage.hpp"
#include "geometry/surface.hpp"
#include "geometry/sine_geometry.hpp"

template<typename value_type, typename fluid_solver, typename time_integrator>
class Swarm : public Surface<Swarm<value_type,fluid_solver,time_integrator> >
{
    public:
        typedef Surface<Swarm<value_type,fluid_solver,time_integrator> > base_type;
        typedef typename base_type::spring_iterator             spring_iterator;
        typedef std::pair<spring_iterator, spring_iterator>     spring_iterator_pair;
        typedef Particle<value_type>                            particle_type;
        typedef SineGeometry<value_type>                        sperm_type;
        typedef std::vector<spring_iterator_pair>               iterator_pair_array;
    private:
        sperm_type m_geometry;
        iterator_pair_array m_tail_iterator_pairs;
        size_t m_num_geometries;

    public:

        Swarm(size_t Mt, size_t Nt, size_t Mh, size_t Nh, int num_sperms = 1)
            : base_type(num_sperms * (Mt * Nt + Mh * (Nh - 1) + 1))
        {
            // Set geometry parameters
            m_geometry.setDimensions(Mt, Nt, Mh, Nh);            
            m_geometry.setWaveSpeed(.001);
            m_geometry.setTailRadius(.05);
            m_geometry.setHeadRadius(.05 * 5);
            m_geometry.setLength(4.0);
            m_geometry.setTailAmplitude(.25);
            m_geometry.setTailPitch(4.1);
            // create a grid where to put the geometries
            std::vector<value_type> mesh2d;
            setGeometryGrid(mesh2d,num_sperms);
            std::vector<size_t> col_ptr, col_idx;
            std::vector<value_type> strenght;
            col_ptr.push_back(0);
            for(size_t i = 0, idx = 0; i < num_sperms; ++i, idx+=3)
            {
                m_geometry.setX0(mesh2d[idx],mesh2d[idx+1],mesh2d[idx+2]);
                m_geometry.init(&this->particles()[i]);
                m_geometry.getConnections(col_ptr,col_idx,i*m_geometry.numParticles());
            }
            getStrengths(col_ptr, col_idx, strenght);
            base_type::setSprings(col_ptr, col_idx, strenght);
        }

        void setIteratorRanges(int num_geometries, int head_offset, int tail_offset)
        {
            spring_iterator s = this->springs_begin(), f;
            
            for (spring_iterator s_end = this->springs_end(); s != s_end; ++s)
                if (s->A()->i % head_offset == 0)
                    break;
                f = s;
            for (spring_iterator s_end = this->springs_end(); f != s_end; ++f)
                if (f->A()->i == hi)
                {
                    do ++f; while (f->A()->i == hi);
                    break;
                }
                m_spring_range = std::make_pair(s, f);
        }
        
        void getStrengths(const std::vector<size_t> &col_ptr, const std::vector<size_t> &col_idx, std::vector<value_type> &strengths)
        {
            strengths.resize(col_idx.size(),10.0);
            for (size_t p = 0, p_end = col_ptr.size() - 1; p < p_end; ++p)
                if (p >= lo && p <= hi)
                    for (size_t i = col_ptr[p], i_end = col_ptr[p + 1]; i < i_end; ++i)
                        strengths[i] = 1.0;
        }
        
        template<typename spring_system_type>
        void set_springs(spring_system_type &spring_system)
        {
            typedef typename spring_system_type::particle_type particle_type;
            particle_type *particles = spring_system.particles();
            size_t *dimensions = m_geometry->get_dimensions();
            size_t tail_offset = dimensions[0] * dimensions[1];
            size_t head_offset = dimensions[2] * (dimensions[3] - 1) + 1;
            size_t offset = tail_offset + head_offset;

            for(size_t i = 0; i < m_num_geometries; ++i)
            {
                size_t idx = i * offset;
                for(size_t p = 0; p < m_head_col_ptr.size() - 1; ++p)
                    for(size_t i = m_head_col_ptr[p], end = m_head_col_ptr[p + 1]; i < end; ++i)
                        if(!spring_system.exist_spring(&particles[p + idx], &particles[m_head_col_idx[i] + idx]))
                            spring_system.add_spring(&particles[p + idx], &particles[m_head_col_idx[i] + idx], .1);

                for(size_t p = 0; p < m_tail_col_ptr.size() - 1; ++p)
                    for(size_t i = m_tail_col_ptr[p], end = m_tail_col_ptr[p + 1]; i < end; ++i)
                        if(!spring_system.exist_spring(&particles[p + idx + head_offset], &particles[m_tail_col_idx[i] + idx + head_offset]))
                            spring_system.add_spring(&particles[p + idx + head_offset], &particles[m_tail_col_idx[i] + idx + head_offset], .1, true);

            }
            spring_system.fluid_solver().initMaps(spring_system);

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
                m_geometry->resetRestingLength(spring, time);
            }
            updateForceGradient(spring_system);
        }


        void setGeometryGrid(std::vector<value_type> &mesh2d, size_t num_geometries)
        {
            value_type dtheta = 1., dalpha = 1.;
            for(size_t i = 0, end = std::sqrt(num_geometries); i < end; ++i)
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
            size_t *dimensions = m_geometry->get_dimensions();
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
            return m_geometry;
        }

};


template<typename _value_type>
struct surface_traits<Swarm<_value_type> >
{
    typedef _value_type value_type;
    typedef SineGeometry<value_type> geometry_type;
};

#endif
