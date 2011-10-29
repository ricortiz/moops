#ifndef SURFACE_HPP
#define SURFACE_HPP

#include <list>
#include <map>

template<typename T> struct surface_traits;

template<typename Derived>
class Surface
{
    protected:
        typedef typename surface_traits<Derived>::geometry_type geometry_type;
        typedef typename surface_traits<Derived>::value_type value_type;
        typedef typename std::map<value_type*,std::pair<size_t,size_t> > grid_type;

    protected:
        grid_type m_grid;

    public:
        Surface() {}
        ~Surface() {}

        inline std::pair<size_t,size_t> &grid(value_type *p)
        {
            return m_grid[p];
        }
        
        inline grid_type &grid()
        {
            return m_grid;
        }
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

        inline geometry_type *geometry() { return derived()->geometry(); }

};


#endif
