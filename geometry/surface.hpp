#ifndef SURFACE_HPP
#define SURFACE_HPP

#include <map>
#include "particle_system/particle.hpp"

template<typename T> struct surface_traits;

template<typename Derived>
class Surface
{
    protected:
        typedef typename surface_traits<Derived>::geometry_type geometry_type;
        typedef typename surface_traits<Derived>::value_type value_type;

    protected:

    public:
        Surface() {}
        ~Surface() {}

        inline Derived *derived()
        {
            return static_cast<Derived*>(this);
        }

        template<typename spring_system_type>
        inline void set_springs(spring_system_type &spring_system)
        {
            derived()->set_springs(spring_system);
        }

        template<typename boundary_type>
        inline void update(boundary_type &boundary, value_type time)
        {
            derived()->update(boundary,time);
        }

};


#endif
