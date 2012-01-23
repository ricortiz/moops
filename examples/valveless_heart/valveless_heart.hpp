#ifndef HEART_PUMP_HPP
#define HEART_PUMP_HPP

#include "geometry/oval_geometry.hpp"
#include "particle_system/elastic_system/elastic_boundary.hpp"
#include <geometry/surface.hpp>


class HeartPump : Surface<HeartPump>
{
    protected:
        typedef ElasticBoundary<typename particle_system_type::particle_type> base_type;
        typedef typename base_type::spring_type spring_type;
        typedef typename base_type::spring_iterator spring_iterator;
	
	typedef typename particle_system_type::value_type value_type;

        typedef OvalGeometry<value_type> oval_type;

    private:
	particle_system_type particle_system;
        fluid_solver_type fluid_solver;
	time_integrator_type time_integrator;
        oval_type &m_geometry;
        std::pair<spring_iterator, spring_iterator> m_spring_range;

    public:

	HeartPump(oval_type &oval_geometry) : particle_system(oval_geometry.numParticles()), fluid_solver(oval_geometry.numParticles()), time_integrator(fluid_solver), m_geometry(oval_geometry)
	{
	  time_integrator.init(particle_system.positions(),particle_system.velocities());
	}
        template<typename particle_type>
        void initSurface(oval_type &oval_geometry, particle_type *particles)
        {
            logger.startTimer("initSurface");
            m_geometry = &oval_geometry;
            m_geometry.init(particles);
            logger.stopTimer("initSurface");
        }

        template<typename particle_type>
        void initVolume(oval_type &oval_geometry, particle_type *particles, size_t num_sub_surfaces = 1)
        {
            oval_geometry.init(particles, num_sub_surfaces);
        }

        void getStrengths(const std::vector<size_t> &col_ptr, const std::vector<size_t> &col_idx, std::vector<value_type> &strengths)
        {
            size_t lo, hi, M, N;
            m_geometry.getForcingRange(lo, hi);
            m_geometry.getDimensions(M, N);
            lo *= M;
            hi *= M;
            for (size_t p = 0; p < col_ptr.size() - 1; ++p)
            {
                if (p > lo && p < hi)
                    for (size_t i = col_ptr[p], end = col_ptr[p + 1]; i < end; ++i)
                        strengths.push_back(1.0);
                else
                    for (size_t i = col_ptr[p], end = col_ptr[p + 1]; i < end; ++i)
                        strengths.push_back(10.0);
            }
        }


        void setSprings()
        {
            size_t lo, hi, M, N;
            m_geometry.getForcingRange(lo, hi);
            m_geometry.getDimensions(M, N);
            std::vector<size_t> col_ptr, col_idx;
            std::vector<value_type> strenght;
            col_ptr.push_back(0);
            m_geometry.getConnections(col_ptr, col_idx);
            getStrengths(col_ptr, col_idx, strenght);

            this->setSprings(col_ptr, col_idx, strenght);
            m_spring_range = getIteratorRange();
        }


        std::pair<spring_iterator, spring_iterator> getIteratorRange()
        {
            size_t lo, hi;
            m_geometry.getForcingRange(lo, hi);
            spring_iterator s, f;
            for (s = this->springs_begin(); s != this->springs_end(); ++s)
                if (s->A()->i == lo || s->B()->i == lo)
                    break;
            f = s;
            for (; f != this->springs_end(); ++f)
                if (f->A()->i == hi)
                {
                    do ++f; while (f->A()->i == hi);
                    break;
                }
            return std::make_pair(s,f);
        }

        inline void updateSprings(value_type time)
        {
            logger.startTimer("updateSprings");
            geometry.setRadiusScaling(time);
            for (spring_iterator s = m_spring_range.first, end = m_spring_range.second; s != end; ++s)
                m_geometry.resetRestingLength(s, time);
            logger.stopTimer("updateSprings");
        }



        inline BaseGeometry<oval_type> &geometry()
        {
            return m_geometry;
        }

};



#endif
