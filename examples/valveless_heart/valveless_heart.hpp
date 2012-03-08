#ifndef HEART_PUMP_HPP
#define HEART_PUMP_HPP

#include "particle_system/surface.hpp"
#include "geometry/oval_geometry.hpp"
#include "particle_system/surface.hpp"

template<typename value_type, typename fluid_solver, typename time_integrator>
class HeartPump : public Surface<HeartPump<value_type, fluid_solver, time_integrator> >
{
    protected:
        typedef Surface<HeartPump<value_type, fluid_solver, time_integrator> > base_type;
        typedef typename base_type::spring_iterator             spring_iterator;
        typedef std::pair<spring_iterator, spring_iterator>     spring_iterator_pair;
        typedef ParticleWrapper<value_type>                     particle_type;
        typedef OvalGeometry<value_type>                        oval_type;

    private:
        oval_type            m_geometry;
        spring_iterator_pair m_spring_range;

    public:
        HeartPump(size_t M, size_t N) : base_type(M*N)
        {
            size_t lo = 10, hi = 50;
            m_geometry.setDimensions(M, N);
            m_geometry.setX0(0, 0, 0);
            m_geometry.setForcingRange(lo, hi);
            m_geometry.setInnerRadius(.05);
            m_geometry.setOuterRadius(.5);
            std::vector<size_t> /*&*/col_ptr/* = this->fluid_solver().col_ptr*/;// sparse matrix (CSR-format) holding
            std::vector<size_t> /*&*/col_idx/* = this->fluid_solver().col_idx*/;// sparse matrix (CSR-format) holding
            std::vector<value_type> strenght;            // interactions between particles
            size_t num_springs = this->particles_size() * 9; // estimate total number of springs
            col_ptr.reserve(this->particles_size() + 1); // Reserve
            col_idx.reserve(num_springs);                // Reserve
            strenght.reserve(num_springs);               // Reserve
            col_ptr.push_back(0);
            m_geometry.init(this->particles());
//             value_type T[3] = {0,0,.15};
//             m_geometry.applyTranslation(this->particles(),T);
            m_geometry.getConnections(col_ptr, col_idx);
            sortConnections(col_ptr, col_idx);
            getStrengths(col_ptr, col_idx, strenght);
            setSprings(col_ptr, col_idx, strenght);
            setIteratorRange(lo, hi);
        }

        inline void computeForces(value_type time)
        {
            m_geometry.setPeristalticRadiusScaling(time);
            for (spring_iterator s = m_spring_range.first, end = m_spring_range.second; s != end; ++s)
                m_geometry.resetRestingLength(s);
            this->clearForces();
            base_type::computeForces();
        }

    private:
        void sortConnections(std::vector<size_t> &col_ptr, std::vector<size_t> &col_idx)
        {
            std::vector<size_t>::iterator begin, end;
            for (size_t i = 0; i < col_ptr.size() - 1; ++i)
            {
                begin = col_idx.begin() + col_ptr[i];
                end = col_idx.begin() + col_ptr[i + 1];
                std::sort(begin, end);
            }
        }

        void getStrengths(const std::vector<size_t> &col_ptr, const std::vector<size_t> &col_idx, std::vector<value_type> &strengths)
        {
            strengths.resize(col_idx.size(), 1.0);            
        }

        void setIteratorRange(size_t lo, size_t hi)
        {
            size_t M, N;
            m_geometry.getDimensions(M, N);
            lo*=M; hi *= M;
            spring_iterator s = this->springs_begin(), s_end = this->springs_end(), f;
            for (; s != s_end; ++s)
                if (s->getAidx() / 3 == lo)
                    break;
            f = s;
            for (; f != s_end; ++f)
                if (f->getAidx() / 3 == hi)
                    break;
            m_spring_range = std::make_pair(s, f);
        }

    public:
        oval_type &geometry()
        {
            return m_geometry;
        }

        template<typename tracers_type>
        void setTracers(tracers_type &tracers, size_t num_rings)
        {
            m_geometry.init(tracers.particles(),num_rings);
//             value_type T[3] = {0,0,.15};
//             m_geometry.applyTranslation(tracers.particles(),T);
        }
};

#include "particle_system/storage/particle_system_storage.hpp"
template<typename _value_type, typename _fluid_solver, typename _time_integrator>
struct Traits<HeartPump<_value_type, _fluid_solver, _time_integrator> >
{
    typedef _value_type value_type;
    typedef _fluid_solver fluid_solver_type;
    typedef _time_integrator time_integrator_type;
    typedef ParticleWrapper<value_type> particle_type;
    typedef ParticleSystemStorage<value_type, particle_type, SURFACE> storage_type;
};


#endif
