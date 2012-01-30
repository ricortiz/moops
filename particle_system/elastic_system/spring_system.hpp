#ifndef SPRING_SYSTEM_HPP
#define SPRING_SYSTEM_HPP

/// @name SpringSystem
/// @section Description
/// @section See also

#include "particle_system/forces/spring.hpp"

template <typename Derived>
class SpringSystem
{
    public:
        typedef typename surface_traits<Derived>::particle_type particle_type;
        typedef typename surface_traits<Derived>::value_type    value_type;
        typedef Spring<particle_type>                           spring_type;
        typedef typename std::list<spring_type>                 spring_container;
        typedef typename spring_container::iterator             spring_iterator;
        typedef typename std::list<spring_type *>               spring_ptr_container;
        typedef typename std::map < particle_type *,
                spring_ptr_container >      spring_lut_type;

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
