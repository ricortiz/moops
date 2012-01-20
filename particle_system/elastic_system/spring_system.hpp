#ifndef SPRING_SYSTEM_HPP
#define SPRING_SYSTEM_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    SpringSystem
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name SpringSystem
/// @section Description
/// @section See also
#include<list>
#include<sstream>
#include<map>
#include<algorithm>
#include "particle_system/forces/spring.hpp"
#include "geometry/surface.hpp"

template <typename surface_type, typename particle_system_type>
class SpringSystem : public particle_system_type
{
    public:

        typedef typename particle_system_type::value_type       value_type;
        typedef typename particle_system_type::particle_type    particle_type;
        typedef Spring<value_type, particle_type>                               spring_type;
        typedef typename std::list<spring_type>                 spring_container;
        typedef typename spring_container::iterator                      spring_iterator;
        typedef typename std::list<spring_type *>                spring_ptr_container;
        typedef typename std::map<particle_type *, spring_ptr_container > spring_lut_type;

    private:
        spring_container        m_springs;      ///< Internal data structure used to store all spring constraints.
        spring_lut_type         m_spring_lut;   ///< Internal datas tructure to record spring connections.
        std::pair<spring_iterator, spring_iterator>    m_spring_update;
        surface_type   *m_surface;

    public:

        SpringSystem() : m_surface(0) {}
        SpringSystem(size_t data_size) : m_surface(0), particle_system_type(data_size) {}
        ~SpringSystem() {}

        inline surface_type *surface()                  { return m_surface; }
        inline surface_type const *surface() const      { return m_surface; }
        spring_lut_type const &springs_map() const      { return m_spring_lut; }
        spring_container const &get_springs() const     { return m_springs; }
        spring_lut_type &springs_map()                  { return m_spring_lut; }
        spring_container &get_springs()                 { return m_springs; }
        spring_iterator springs_begin()                 { return m_springs.begin(); }
        spring_iterator springs_end()                   { return m_springs.end(); }
        size_t springs_size()                           { return m_springs.size(); }
        void set_surface(surface_type &surface)         { m_surface = &surface; }

        inline void computeForces()
        {
            this->clear_forces();
            for(spring_iterator s = m_springs.begin(), end = m_springs.end(); s != end; ++s)
                s->apply();
        }

        inline void computeForces(const value_type *x)
        {
	  this->clear_forces();
            for(spring_iterator s = m_springs.begin(), end = m_springs.end(); s != end; ++s)
            {
                size_t i = s->getAidx();
                size_t j = s->getBidx();
                const value_type *x1 = &x[i];
                const value_type *x2 = &x[j];
                s->apply(x1, x2);
            }
        }

        inline void computeForces(const value_type *x, value_type *f)
        {
	  std::fill(f,f+this->data_size(),0.0);
            for(spring_iterator s = m_springs.begin(), end = m_springs.end(); s != end; ++s)
            {
                size_t i = s->getAidx();
                size_t j = s->getBidx();
                value_type *x1 = &x[i];
                value_type *x2 = &x[j];
                value_type *f1 = &f[i];
                value_type *f2 = &f[j];
                s->apply(x1, x2, f1, f2);
            }
        }

        inline void updateForces(value_type time, const value_type *x, value_type *f)
        {
            if(m_surface)
                m_surface->updateSprings(*this, time);
            computeForces(x, f);
        }

        inline void updateForces(value_type time, const value_type *x)
        {	  
            if(m_surface)
                m_surface->updateSprings(*this, time);
            computeForces(x);
        }
        
        inline void updateForces(value_type time)
        {
            if(m_surface)
                m_surface->updateSprings(*this, time);
            computeForces();
        }
        
        void clear()
        {
            m_spring_lut.clear();
            m_springs.clear();
        }

        inline spring_iterator addSpring(particle_type *A, particle_type *B, value_type k = 1)
        {
            m_springs.push_back(spring_type());
            spring_iterator s = m_springs.end();
            --s;

            s->init(A, B);
            s->stiffness() = k;
            m_spring_lut[&(*A)].push_back(&(*s));
            m_spring_lut[&(*B)].push_back(&(*s));

            return s;
        }

        inline bool existSpring(particle_type *A, particle_type *B)
        {
            typedef typename spring_ptr_container::iterator spring_ptr_iterator;
            if(A == B)
                return true;
            spring_ptr_container springsA = m_spring_lut[A];
            spring_ptr_container springsB = m_spring_lut[B];
            if(springsA.empty())
                return false;
            if(springsB.empty())
                return false;
            {
                spring_ptr_iterator  s   = springsA.begin();
                spring_ptr_iterator  end = springsA.end();
                for(; s != end; ++s)
                {
                    if(((*s)->A() == A && (*s)->B() == B) || ((*s)->B() == A && (*s)->A() == B))
                        return true;
                }
            }
            {
                spring_ptr_iterator  s   = springsB.begin();
                spring_ptr_iterator  end = springsB.end();
                for(; s != end; ++s)
                {
                    if(((*s)->A() == A && (*s)->B() == B) || ((*s)->B() == A && (*s)->A() == B))
                        return true;
                }
            }
            return false;
        }

        template<typename int_vector_type, typename real_vector_type>
        inline void setSprings(int_vector_type &col_ptr, int_vector_type &col_idx, real_vector_type &strength, size_t min = 0, size_t max = 0)
        {
            logger.startTimer("setSprings");
            particle_type *particles = this->particles();
            for(size_t p = 0; p < col_ptr.size() - 1; ++p)
                for(size_t i = col_ptr[p], end = col_ptr[p + 1]; i < end; ++i)
                    if(!existSpring(&particles[p], &particles[col_idx[i]]))
                    {
                        spring_iterator s = addSpring(&particles[p], &particles[col_idx[i]], strength[i]);
                        s->getAidx() = 3 * p;
                        s->getBidx() = 3 * col_idx[i];
                    }

            spring_iterator p;
            for(p = m_springs.begin(); p != m_springs.end(); ++p)
                if(p->A()->i == min || p->B()->i == min)
                {
                    m_spring_update.first = p;
                    break;
                }
            for(; p != m_springs.end(); ++p)
                if(p->A()->i == max)
                {
                    do ++p; while(p->A()->i == max);
                    break;
                }

            m_spring_update.second = p;

            logger.stopTimer("setSprings");
        }

        template<typename geometry_type>
        inline void updateSprings(geometry_type &geometry, value_type time)
        {
            logger.startTimer("updateSprings");
            geometry.setRadiusScaling(time);
// 	    std::cout << "springs = [";
// 	    for(spring_iterator s = m_spring_update.first, end = m_spring_update.second; s != end; ++s)
// 	      std::cout << s->getAidx()/3 << "," << s->getBidx()/3 << ";";
// 	    std::cout << "];" << std::endl;
// 	    std::stringstream s1,s2;
// 	    s1 << "lengths1 = [";
// 	    s2 << "lengths2 = [";
            for(spring_iterator s = m_spring_update.first, end = m_spring_update.second; s != end; ++s)
	    {
// 	      s1 << s->resting_length() << " ";
                geometry.resetRestingLength(s, time);
// 	      s2 << s->resting_length() << " ";
	    }
// 	    s1 << "];" << std::endl;
// 	    s2 << "];" << std::endl;
// 	    std::cout << s1.str();
// 	    std::cout << s2.str();
	    
            logger.stopTimer("updateSprings");
        }

};
#endif
