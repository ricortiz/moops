#ifndef HEART_PUMP_HPP
#define HEART_PUMP_HPP

#include "geometry/oval_geometry.hpp"
#include "particle_system/storage/particle_system_storage.hpp"
#include "geometry/surface.hpp"

template<typename value_type, typename fluid_solver, typename time_integrator>
class HeartPump : public Surface<HeartPump<value_type, fluid_solver, time_integrator> >
{
    protected:
        typedef Surface<HeartPump<value_type, fluid_solver, time_integrator> > base_type;
        typedef typename base_type::spring_iterator spring_iterator;
        typedef OvalGeometry<value_type> oval_type;
        typedef Particle<value_type> particle_type;
        typedef std::pair<spring_iterator, spring_iterator> spring_iterator_pair;

    private:
        oval_type            m_geometry;
        spring_iterator_pair m_spring_range;

    public:
        HeartPump(size_t M, size_t N) : base_type(M*N)
        {
            size_t lo = 10, hi = 50;
            m_geometry.setDimensions(M,N);
            m_geometry.setX0(0,0,0);
            m_geometry.setForcingRange(lo,hi);
            m_geometry.setWaveSpeed(.001);
            m_geometry.init(this->particles());
            setSprings(lo*M,hi*M);
            setIteratorRange(lo,hi);
        }

        inline void computeForces(value_type time)
        {
            m_geometry.setRadiusScaling(time);
            for (spring_iterator s = m_spring_range.first, end = m_spring_range.second; s != end; ++s)
                m_geometry.resetRestingLength(s);
            base_type::computeForces();
        }

    private:
        void setSprings(size_t plo, size_t phi)
        {
            std::vector<size_t> col_ptr, col_idx;
            std::vector<value_type> strenght;
            col_ptr.push_back(0);
            m_geometry.getConnections(col_ptr, col_idx);
            getStrengths(col_ptr, col_idx, strenght, plo, phi);
            
            base_type::setSprings(col_ptr, col_idx, strenght);
        }
        
        void initVolume(oval_type &oval_geometry, particle_type *particles, size_t num_sub_surfaces = 1)
        {
            oval_geometry.init(particles, num_sub_surfaces);
        }

        void getStrengths(const std::vector<size_t> &col_ptr, const std::vector<size_t> &col_idx, std::vector<value_type> &strengths, size_t lo, size_t hi)
        {
            strengths.resize(col_idx.size(),10.0);
            for (size_t p = 0, p_end = col_ptr.size() - 1; p < p_end; ++p)
                if (p >= lo && p <= hi)
                    for (size_t i = col_ptr[p], i_end = col_ptr[p + 1]; i < i_end; ++i)
                        strengths[i] = 1.0;            
        }

        void setIteratorRange(size_t lo, size_t hi)
        {
            spring_iterator s = this->springs_begin(), f;
            for (spring_iterator s_end = this->springs_end(); s != s_end; ++s)
                if (s->A()->i == lo || s->B()->i == lo)
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

};

template<typename _value_type, typename _fluid_solver, typename _time_integrator>
struct surface_traits<HeartPump<_value_type, _fluid_solver, _time_integrator> >
{
    typedef _value_type value_type;
    typedef _fluid_solver fluid_solver_type;
    typedef _time_integrator time_integrator_type;
    typedef Particle<value_type> particle_type;
    typedef ParticleSystemStorage<value_type, particle_type, SURFACE> storage_type;
};


#endif
