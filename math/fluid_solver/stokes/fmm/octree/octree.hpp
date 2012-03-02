#ifndef OCTREE_HPP
#define OCTREE_HPP
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
#include<cmath>
#include <omp.h>
#include "math/fluid_solver/stokes/fmm/octree/octree_storage.hpp"
#include "math/fluid_solver/stokes/fmm/octree/box.hpp"

/**
 * @brief The Octree class. This class creates the octree data structure.
 *
 * @param value_type Real type, eg. float or double
 * @param box_type Box type of each node in the tree
 * @param max_particles This is the maximum number of particles allowed in each box
 * @param array_type This is the internal storage type.  It is used to allocate and store the boxes.
 *
 **/
template<typename value_type, size_t max_particles, typename particle_type, typename box_type = Box<value_type,particle_type>, typename array_type = std::list<box_type> >
class Octree
{
    protected:
        typedef typename std::list<box_type*>::iterator            box_iterator;
        typedef typename std::list<box_type*>::const_iterator      const_box_iterator;

    private:
        octree_storage<box_type, array_type> m_storage;
        size_t m_depth;
        particle_type *m_particles;

    public:
        inline Octree(value_type center[], value_type extent[]) : m_depth(0), m_storage()
        {
            root()->x() = center[0];
            root()->y() = center[1];
            root()->z() = center[2];
            root()->u() = extent[0];
            root()->v() = extent[1];
            root()->w() = extent[2];
        }
        
        inline Octree(value_type *positions, value_type *velocities, value_type *forces, size_t num_particles)
        {
            m_particles = new particle_type[num_particles];
            link_particles(positions,velocities,forces,m_particles,num_particles);
            get_dimensions(m_particles,num_particles,root()->coords(),root()->extent());
            set_root_particles(m_particles,num_particles);
        }
        
        inline Octree(particle_type *particles, size_t num_particles) : m_particles(particles)
        {
            get_dimensions(m_particles,num_particles,root()->coords(),root()->extent());
            set_root_particles(m_particles,num_particles);
        }
        ~Octree() { }

        inline box_type *root()                { return &m_storage.root(); }
        inline box_type const *root() const    { return &m_storage.root(); }
        inline const array_type &boxes() const { return m_storage.boxes(); }
        inline array_type &boxes()             { return m_storage.boxes(); }
        inline size_t &depth()                 { return m_depth; }
        inline size_t const &depth() const     { return m_depth; }
        
        /**
         * @brief Create a subtree at provided node box.
         *
         * @param box Node where to start a subtree
         * 
         **/
        inline void build_tree(box_type *box)
        {
            if (box->particles().size() > max_particles)
            {
                box->create_children(boxes());
                size_t depth = boxes().back().level()+1;
                if(depth > m_depth)
                    m_depth = depth;
            }
            for (box_iterator b = box->children().begin(), end = box->children().end(); b != end; ++b)
            {
#pragma omp task firstprivate(b)
                {
                    build_tree(*b);
                }
            }
#pragma omp taskwait
        }

        /**
         * @brief Build octree from root node.
         *
         * 
         **/
        inline void build_tree()
        {
            build_tree(root());
        }

        
        /**
         * @brief Copy particle array into the root node of the tree.  It actually only copy a pointer to the particles
         *
         * @param particles particle array
         * @param num_particles total # of particles
         * 
         **/
        inline void set_root_particles(particle_type particles[], size_t num_particles)
        {
            std::list<particle_type*> &particle_array = root()->particles();
            for (size_t i = 0; i < num_particles; ++i)
                particle_array.push_back(&particles[i]);
        }
        
        inline void link_particles(value_type *p, value_type *v, value_type *f, particle_type *particles, size_t num_particles)
        {
            for (size_t i = 0, idx = 0; i < num_particles; ++i, idx+=3)
            {
                particles[i].position = &p[idx];
                particles[i].velocity = &v[idx];
                particles[i].force = &f[idx];
                
            }
        }
        
        inline void get_dimensions(particle_type particles[], size_t num_particles, value_type center[], value_type extent[])
        {
            value_type min[3] = {1e100,1e100,1e100};
            value_type max[3] = {-1e100,-1e100,-1e100};
            center[0] = center[1] = center[2] = 0;
            extent[0] = extent[1] = extent[2] = 0;
            for (size_t i = 0; i < num_particles; i++)
            {
                
                for (int k = 0 ; k < 3; ++k)
                {
                    if (particles[i].position[k] > max[k])
                        max[k] = particles[i].position[k];
                    if (particles[i].position[k] < min[k])
                        min[k] = particles[i].position[k];
                    center[k] += particles[i].position[k];
                }
            }
            
            for (int i = 0; i < 3; ++i)
            {
                center[i] /= num_particles;
                max[i] = std::abs(max[i]-center[i])+.001;
                min[i] = std::abs(min[i]-center[i])+.001;
                extent[i] = std::max(max[i],min[i]);
            }
        }
        
        /**
         * @brief Build tree and initialize lists
         * 
         **/
        inline void init()
        {
#pragma omp parallel
            {
#pragma omp single
                {
                    build_tree();
                    root()->buildColleagues();
                    root()->buildList1();
                    root()->buildList2();
                    root()->buildList3();
                    root()->buildList4();
                }
            }
#pragma omp barrier
        }

        
        /**
         * @brief Overloaded output operator
         *
         * @param o output stream
         * @param tree tree to output
         * @return output_type& output stream
         **/
        template<typename output_type>
        friend output_type & operator<< (output_type & o, Octree const &tree)
        {
            typedef typename array_type::const_iterator iterator;
            std::cout << "Tree depth: " << tree.m_depth << std::endl;
            std::cout << "Num boxes: " << tree.boxes().size() << std::endl;
            for (iterator b = tree.boxes().begin(), end = tree.boxes().end(); b != end; ++b)
                o << *b << std::endl;
            return o;
        }

};


#endif


