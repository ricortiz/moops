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
#include<map>
#include "particle_system/forces/spring.hpp"

template <typename Derived>
class SpringSystem
{
    public:
        surface_traits<Derived>::particle_type particle_type;
        surface_traits<Derived>::value_type value_type;
        typedef Spring<particle_type>                               spring_type;
        typedef typename std::list<spring_type>                 spring_container;
        typedef typename spring_container::iterator                      spring_iterator;
        typedef typename std::list<spring_type *>                spring_ptr_container;
        typedef typename std::map<particle_type *, spring_ptr_container > spring_lut_type;

    private:
        spring_container        m_springs;      ///< Internal data structure used to store all spring constraints.
        spring_lut_type         m_spring_lut;   ///< Internal datas tructure to record spring connections.

    public:

        SpringSystem() {}
        ~SpringSystem() {}

        spring_lut_type const &springs_map() const      { return m_spring_lut; }
        spring_container const &get_springs() const     { return m_springs; }
        spring_lut_type &springs_map()                  { return m_spring_lut; }
        spring_container &get_springs()                 { return m_springs; }
        spring_iterator springs_begin()                 { return m_springs.begin(); }
        spring_iterator springs_end()                   { return m_springs.end(); }
        size_t springs_size()                           { return m_springs.size(); }

        inline void computeForces()
        {
            for (spring_iterator s = m_springs.begin(), end = m_springs.end(); s != end; ++s)
                s->apply();
        }

        inline void computeForces(const value_type *x)
        {
            for (spring_iterator s = m_springs.begin(), end = m_springs.end(); s != end; ++s)
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
            for (spring_iterator s = m_springs.begin(), end = m_springs.end(); s != end; ++s)
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
            if (A == B)
                return true;
            spring_ptr_container springsA = m_spring_lut[A];
            spring_ptr_container springsB = m_spring_lut[B];
            if (springsA.empty())
                return false;
            if (springsB.empty())
                return false;
            {
                spring_ptr_iterator  s   = springsA.begin();
                spring_ptr_iterator  end = springsA.end();
                for (; s != end; ++s)
                {
                    if (((*s)->A() == A && (*s)->B() == B) || ((*s)->B() == A && (*s)->A() == B))
                        return true;
                }
            }
            {
                spring_ptr_iterator  s   = springsB.begin();
                spring_ptr_iterator  end = springsB.end();
                for (; s != end; ++s)
                {
                    if (((*s)->A() == A && (*s)->B() == B) || ((*s)->B() == A && (*s)->A() == B))
                        return true;
                }
            }
            return false;
        }



};
#endif
