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

    private:
        oval_type *m_geometry;
        
    public:

        template<typename particle_type>
        void initSurface(oval_type &oval_geometry, particle_type *particles)
        {
            logger.startTimer("initSurface");
            m_geometry = &oval_geometry;
            m_geometry->init(particles);
            logger.stopTimer("initSurface");
        }
        
        template<typename particle_type>
        void initVolume(BaseGeometry<oval_type> &oval_geometry, particle_type *particles, size_t num_sub_surfaces = 1)
        {
            oval_geometry.init(particles, num_sub_surfaces);
        }

        void setStrengths(std::vector<size_t> &col_ptr, std::vector<size_t> &col_idx, std::vector<value_type> &strengths)
        {
            size_t lo, hi, M, N;
            m_geometry->getForcingRange(lo,hi);
            m_geometry->getDimensions(M,N);
            lo *= M;
            hi *= M;
            for(size_t p = 0; p < col_ptr.size() - 1; ++p)
            {
                if(p > lo && p < hi)
                    for(size_t i = col_ptr[p], end = col_ptr[p + 1]; i < end; ++i)
                        strengths.push_back(1.0);
                else
                    for(size_t i = col_ptr[p], end = col_ptr[p + 1]; i < end; ++i)
                        strengths.push_back(10.0);
            }
        }

        template<typename spring_system_type>
        void setSprings(spring_system_type &spring_system)
        {
            size_t lo, hi, M, N;
            m_geometry->getForcingRange(lo,hi);
            m_geometry->getDimensions(M,N);
            std::vector<size_t> col_ptr, col_idx;
            std::vector<value_type> strenght;
            col_ptr.push_back(0);
            m_geometry->getConnections(col_ptr, col_idx);
            setStrengths(col_ptr,col_idx,strenght);
            
            spring_system.setSprings(col_ptr,col_idx,strenght,lo,hi);
        }

        template<typename spring_system_type>
        void updateSprings(spring_system_type &spring_system, value_type time)
        {
            spring_system.updateSprings(*m_geometry,time);
        }

        inline BaseGeometry<oval_type> &geometry()
        {
            return *m_geometry;
        }

};


template<typename _value_type>
struct surface_traits<HeartPump<_value_type> >
{
    typedef _value_type value_type;
    typedef OvalGeometry<value_type> geometry_type;
};

#endif
