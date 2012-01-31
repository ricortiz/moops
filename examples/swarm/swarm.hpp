#ifndef SWARM_HPP
#define SWARM_HPP

#include "particle_system/particle.hpp"
#include "particle_system/surface.hpp"
#include "geometry/sine_geometry.hpp"

template<typename value_type, typename fluid_solver, typename time_integrator>
class Swarm : public Surface<Swarm<value_type, fluid_solver, time_integrator> >
{
    public:
        typedef Surface<Swarm<value_type, fluid_solver, time_integrator> >      base_type;
        typedef typename base_type::spring_iterator                             spring_iterator;
        typedef std::pair<spring_iterator, spring_iterator>                     spring_iterator_pair;
        typedef Particle<value_type>                                            particle_type;
        typedef SineGeometry<value_type>                                        sperm_type;
        typedef std::vector<spring_iterator_pair>                               iterator_pair_array;

    private:
        sperm_type              m_geometry;             //< Contains geometrical props of sperms
        iterator_pair_array     m_tail_iterator_pairs;  //< Stores spring iterator range for tail
        int                     m_num_geometries;

    public:

        Swarm(size_t Mt, size_t Nt, size_t Mh, size_t Nh, int num_sperms = 1)
            : base_type(num_sperms * (Mt * Nt + Mh * (Nh - 1) + 1)), m_num_geometries(num_sperms)
        {
            // Set geometry parameters
            m_geometry.setDimensions(Mt, Nt, Mh, Nh);    // tail dims: MtxNt; head dims: MhxNh
            m_geometry.setWaveSpeed(.001);               // speed of wave passed to tail
            m_geometry.setTailRadius(.05);               // radius of tail tube
            m_geometry.setHeadRadius(.05 * 5);           // radius of head
            m_geometry.setLength(4.0);                   // length of the tail
            m_geometry.setTailAmplitude(.25);            // initial amplitude of tail
            m_geometry.setTailPitch(4.1);                // pitch of tail
            std::vector<value_type> mesh2d;              // coords of each geometry
            size_t x = std::sqrt(num_sperms) , y = x;
            setGeometryGrid(mesh2d, x, y);         // create a grid to put the geometries on
            std::vector<size_t> col_ptr, col_idx;        // sparse matrix (CSR-format) holding
            std::vector<value_type> strenght;            // interactions between particles
            size_t num_springs = this->particles_size() * 9; // estimate total number of springs
            col_ptr.reserve(this->particles_size() + 1); // Reserve
            col_idx.reserve(num_springs);                // Reserve
            strenght.reserve(num_springs);               // Reserve
            col_ptr.push_back(0);
            m_geometry.init(&this->particles()[0]);
            m_geometry.getConnections(col_ptr, col_idx);
            for (int i = 1, idx = 3; i < num_sperms; ++i, idx += 3)
            {
                particle_type *p_init = &this->particles()[0];
                particle_type *p = &this->particles()[i * m_geometry.numParticles()];
                for (size_t j = 0; j < m_geometry.numParticles(); ++j)
                {
                    p[j].position[0] = p_init[j].position[0] + mesh2d[idx];
                    p[j].position[1] = p_init[j].position[1] + mesh2d[idx + 1];
                    p[j].position[2] = p_init[j].position[2] + mesh2d[idx + 2];
                    p[j].i = p_init[j].i;
                    p[j].j = p_init[j].j;
                }
                m_geometry.getConnections(col_ptr, col_idx, i * m_geometry.numParticles());
            }
            getStrengths(col_ptr, col_idx, strenght);
            base_type::setSprings(col_ptr, col_idx, strenght);
            
            setIteratorRanges(Mt * Nt, Mh * (Nh - 1) + 1);
	    
	    this->fluid_solver().initMaps(*this);
        }

        inline void computeForces(value_type time)
        {
            for(size_t i = 0; i < m_tail_iterator_pairs.size(); ++i)
            {
                for (spring_iterator s = m_tail_iterator_pairs[i].first, end = m_tail_iterator_pairs[i].second; s != end; ++s)
                    m_geometry.resetRestingLength(s, time);
            }
            this->clear_forces();
            base_type::computeForces();
            updateForceGradient();
        }

    private:
        void setGeometryGrid(std::vector<value_type> &mesh2d, int x, int y)
        {
            value_type dtheta = 1., dalpha = 1.;
            for (int i = 0; i < x; ++i)
                for (int j = 0; j < y; ++j)
                {
                    mesh2d.push_back(i * dtheta);
                    mesh2d.push_back(0.0);
                    mesh2d.push_back(j * dalpha);
                }
        }

        void setIteratorRanges(int tail_offset, int head_offset)
        {
            spring_iterator s = this->springs_begin(), s_end = this->springs_end(), f;
            size_t offset = 0;
            while (s != s_end)
            {
                offset += head_offset;
                for (; s != s_end; ++s)
                    if (s->getAidx() / 3 == offset)
                        break;
                offset += tail_offset;
                f = s;
                for (; f != s_end; ++f)
                    if (f->getAidx() / 3 == offset)
                        break;
                m_tail_iterator_pairs.push_back(std::make_pair(s, f));
                s = f;
            }

        }

        void getStrengths(const std::vector<size_t> &, const std::vector<size_t> &col_idx, std::vector<value_type> &strengths)
        {
            strengths.resize(col_idx.size(), 5.0);
        }

        void updateForceGradient()
        {
            size_t Mt, Nt, Mh, Nh;
            m_geometry.getDimensions(Mt, Nt, Mh, Nh);
            size_t tail_offset = Mt * Nt;
            size_t head_offset = Mh * (Nh - 1) + 1;
            size_t offset = tail_offset + head_offset;

            for(int i = 0; i < m_num_geometries; ++i)
            {
                value_type gradient[3] = {0};
                particle_type *particles = this->particles() + i * offset;
                for (size_t j = 0; j < head_offset; ++j)
                {
                    gradient[0] += particles[j].force[0];
                    gradient[1] += particles[j].force[1];
                    gradient[2] += particles[j].force[2];
                }
                gradient[0] /= head_offset;
                gradient[1] /= head_offset;
                gradient[2] /= head_offset;
//                 value_type norm = std::sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1] + gradient[2] * gradient[2]);

                value_type scale = (1.0/std::sqrt(particles[0].position[0]*particles[0].position[0]*particles[0].position[1]*particles[0].position[1]*particles[0].position[2]*particles[0].position[2]));

                if (scale > 1)
                    return;
                
                gradient[0] = -particles[0].position[0]*scale;
                gradient[1] = -particles[0].position[1]*scale;
                gradient[2] = -particles[0].position[2]*scale;
                particles[0].force[0] += gradient[0];
                particles[0].force[1] += gradient[1];
                particles[0].force[2] += gradient[2];
            }

        }

    public:
        sperm_type &geometry()
        {
            return m_geometry;
        }

        template<typename out_stream>
        void print_springs(out_stream &out = std::cout)
        {
            out << "springs = [";
            spring_iterator s = this->springs_begin(), end = this->springs_end();
            for(; s != end; ++s)
                out << s->getAidx() / 3 + 1 << "," << s->getBidx() / 3 + 1 << ";";
            out << "];" << std::endl;
        }
        
        template<typename out_stream>
        void print_forces(out_stream &out = std::cout)
        {
            value_type *p = this->forces();
            out << "f = [";
            for(size_t i = 0, idx = 0; i < this->particles_size(); ++i, idx += 3)
                out << p[idx] << "," << p[idx + 1] << "," << p[idx + 2] << ";";
            out << "];" << std::endl;
        }

        template<typename out_stream>
        void print_velocities(out_stream &out = std::cout)
        {
            value_type *p = this->velocities();
            out << "v = [";
            for(size_t i = 0, idx = 0; i < this->particles_size(); ++i, idx += 3)
                out << p[idx] << "," << p[idx + 1] << "," << p[idx + 2] << ";";
            out << "];" << std::endl;
        }
        
        template<typename out_stream>
        void print_positions(out_stream &out = std::cout)
        {
            value_type *p = this->positions();
            out << "p = [";
            for(size_t i = 0, idx = 0; i < this->particles_size(); ++i, idx += 3)
                out << p[idx] << "," << p[idx + 1] << "," << p[idx + 2] << ";";
            out << "];" << std::endl;
        }

        template<typename out_stream>
        void print_tail_springs(out_stream &out = std::cout)
        {
            out << "tail_springs = [";
            for(size_t i = 0; i < m_tail_iterator_pairs.size(); ++i)
            {
                spring_iterator s = m_tail_iterator_pairs[i].first, end = m_tail_iterator_pairs[i].second;
                for(; s != end; ++s)
                    out << s->getAidx() / 3 + 1 << "," << s->getBidx() / 3 + 1 << ";";
            }
            out << "];" << std::endl;
        }

        template<typename out_stream>
        void print_spring_lengths(out_stream &out = std::cout)
        {
            out << "spring_lenghts = [";
            for(size_t i = 0; i < m_tail_iterator_pairs.size(); ++i)
            {
                spring_iterator s = m_tail_iterator_pairs[i].first, end = m_tail_iterator_pairs[i].second;
                for(; s != end; ++s)
                    out << s->resting_length() << ";";
            }
            out << "];" << std::endl;
        }

};

#include "particle_system/storage/particle_system_storage.hpp"
template<typename _value_type, typename _fluid_solver, typename _time_integrator>
struct Traits<Swarm<_value_type, _fluid_solver, _time_integrator> >
{
    typedef _value_type value_type;
    typedef _fluid_solver fluid_solver_type;
    typedef _time_integrator time_integrator_type;
    typedef Particle<value_type>                                      particle_type;
    typedef ParticleSystemStorage<value_type, particle_type, SURFACE> storage_type;
};

#endif
