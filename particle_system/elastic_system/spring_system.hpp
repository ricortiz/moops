#ifndef SPRING_SYSTEM_HPP
#define SPRING_SYSTEM_HPP
/****************************************************************************
** MOOPS -- Modular Object Oriented Particle Simulator
** Copyright (C) 2011-2012  Ricardo Ortiz <ortiz@unc.edu>
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/

#include "particle_system/forces/spring.hpp"

template <typename Derived>
class SpringSystem
{
    public:
        typedef typename Traits<Derived>::particle_type particle_type;
        typedef typename Traits<Derived>::value_type    value_type;
        typedef Spring<particle_type>                           spring_type;
        typedef typename std::list<spring_type>                 spring_container;
        typedef typename spring_container::iterator             spring_iterator;
        typedef typename std::list<spring_type*>              spring_ptr_container;
        typedef typename std::map <particle_type*, spring_ptr_container>      spring_lut_type;

    private:
        spring_container        m_springs;      ///< Internal data structure used to store all spring constraints.
        spring_lut_type         m_spring_lut;   ///< Internal datas tructure to record spring connections.

    public:
        inline Derived &derived()
        {
            return *static_cast<Derived*>(this);
        }

        spring_iterator springs_begin()                 { return m_springs.begin(); }
        spring_iterator springs_end()                   { return m_springs.end(); }
        const spring_container &springs() const { return m_springs; }
        spring_container &springs() { return m_springs; }

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

        template<typename out_stream>
        void printSprings(out_stream &out = std::cout)
        {
            printSprings(out, m_springs.begin(), m_springs.end());
        }

        template<typename out_stream>
        void printRestingLengths(out_stream &out = std::cout)
        {
            printRestingLengths(out, m_springs.begin(), m_springs.end());
        }
        
        template<typename out_stream>
        void printSprings(out_stream &out, spring_iterator start, spring_iterator end)
        {
            out << "springs = [";
            spring_iterator s = start;
            for (; s != end; ++s)
                out << s->getAidx() / 3 + 1 << "," << s->getBidx() / 3 + 1 << ";";
            out << "];" << std::endl;
        }

        template<typename out_stream>
        void printRestingLengths(out_stream &out, spring_iterator start, spring_iterator end)
        {
            out << "spring_lenghts = [";
            spring_iterator s = start;
            for (; s != end; ++s)
                out << s->resting_length() << ";";
            out << "];" << std::endl;
        }

        template<typename out_stream>
        friend out_stream &operator<<(out_stream &out, SpringSystem<Derived> &system)
        {
            system.printSprings(out);
            return out;
        }

};
#endif
