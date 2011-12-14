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
#include "geometry/surface.hpp"

template <typename surface_type, typename particle_system_type>
class SpringSystem : public particle_system_type
{
    public:

        typedef typename particle_system_type::value_type       value_type;
        typedef typename particle_system_type::particle_type    particle_type;
        typedef typename particle_system_type::particle_integrator_type    particle_integrator_type;
        typedef Spring<value_type, particle_type>                               spring_type;
        typedef typename std::list<spring_type>                 spring_container;
        typedef typename spring_container::iterator                      spring_iterator;
        typedef typename std::list<spring_type *>                spring_ptr_container;
        typedef typename std::map<particle_type *, spring_ptr_container > spring_lut_type;        
        std::map<value_type*,std::vector<value_type*> > connections;
        
    private:
        spring_container    m_springs;      ///< Internal data structure used to store all spring constraints.
        spring_lut_type     m_spring_lut;   ///< Internal datas tructure to record spring connections.
        Surface<surface_type>        *m_surface;

    public:

        SpringSystem() : m_surface(0) {}
        ~SpringSystem() {}

        inline Surface<surface_type> *surface() { return m_surface; }
        inline Surface<surface_type> const *surface() const { return m_surface; }

        spring_lut_type const &springs_map() const { return m_spring_lut; }
        spring_container const &get_springs() const { return m_springs; }
        spring_lut_type &springs_map() { return m_spring_lut; }
        spring_container &get_springs() { return m_springs; }
        spring_iterator springs_begin() { return m_springs.begin(); }
        spring_iterator springs_end() { return m_springs.end(); }
        void set_surface(Surface<surface_type> &surface) { m_surface = &surface; }
        size_t springs_size() { return m_springs.size(); }

        inline void compute_forces()
        {
            this->clear_forces();
            for(spring_iterator f = m_springs.begin(), end = m_springs.end(); f != end; ++f)
                f->apply();
        }

        inline void update_forces(value_type time)
        {
            if(m_surface)
                m_surface->update(*this, time);
            compute_forces();
        }

        void clear()
        {
            m_spring_lut.clear();
            m_springs.clear();
        }

        void set_resting_lengths(value_type resting_lenghts[])
        {
            size_t i = 0;
            for(spring_iterator s = m_springs.begin(), end = m_springs.end(); s != end; ++s)
                s->resting_length() = resting_lenghts[i++];
        }

        inline bool add_spring(particle_type *A, particle_type *B, value_type k = 1, bool add = true)
        {
            if(A == B)
                return false;
            m_springs.push_back(spring_type());
            spring_type *s = &m_springs.back();
            s->init(A, B);
            s->stiffness() = k;
            connections[A->position].push_back(B->position);
            if(add)
            {
                m_spring_lut[&(*A)].push_back(s);
                m_spring_lut[&(*B)].push_back(s);
            }
            return true;
        }

        inline bool exist_spring(particle_type *A, particle_type *B)
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

};
#endif
