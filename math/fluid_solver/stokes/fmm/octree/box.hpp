#ifndef BOX_HPP
#define BOX_HPP

#include<cassert>
#include<list>

#include "math/fluid_solver/stokes/fmm/utils/bits.hpp"
#include "math/fluid_solver/stokes/fmm/utils/meta.hpp"

/**
 * @brief The Box class.  This class represent a node in a octree.
 *
 * You should not use this class by itself.  It is used passes to the \ref Octree class.
 *
 * @param value_type Real type, eg. float or double
 * @param particle_type The type of the particle containing the location
 * @param dim The dimension of the box, defaults to 3D
 * @param num_children the number of children for the box,  this parameter depends on the dimension.
 *
 * \code
 * struct particle_type { double position[3]; };
 *
 * typedef Box<double,particle_type> box_type;
 *
 * \endcode
 **/
template<typename value_type, typename particle_type, int dim = 3, int num_children = 1 << dim>
class Box
{
    public:
        typedef Box<value_type,particle_type,dim,num_children>     box_type;
        typedef typename std::list<box_type*>::iterator            box_iterator;
        typedef typename std::list<box_type*>::const_iterator      const_box_iterator;
        typedef typename std::list<particle_type*>::iterator       particle_iterator;
        typedef typename std::list<particle_type*>::const_iterator const_particle_iterator;

    private:
        enum
        {
            m_dim = dim,
            m_num_children = num_children
        };
        box_type *m_parent;
        size_t m_level;
        size_t m_id;
        size_t m_idx;
        value_type m_coords[m_dim];
        value_type m_extent[m_dim];
        std::list<particle_type*> m_particles;
        std::list<box_type*> m_children;
        std::list<box_type*> m_colleagues;
        std::list<box_type*> m_lists[4];

    public:
        Box() : m_parent(0), m_level(0), m_id(0), m_idx(0)
        {
            m_coords[0] = m_coords[1] = m_coords[2] = value_type(0);
            m_extent[0] = m_extent[1] = m_extent[2] = value_type(1);
        }
        ~Box() {}

        /// Data access routines
        inline box_type *parent()                       { return m_parent; }
        inline box_type const *parent() const           { return m_parent; }
        inline void set_parent(box_type *box)           { m_parent = box; }
        inline size_t const &level() const              { return m_level; }
        inline size_t &level()                          { return m_level; }
        inline size_t const &id() const                 { return m_id; }
        inline size_t &id()                             { return m_id; }
        inline size_t const &idx() const                 { return m_idx; }
        inline size_t &idx()                             { return m_idx; }
        inline value_type *coords()                   { return m_coords; }
        inline value_type const *coords() const       { return m_coords; }
        inline value_type &coords(int i)             { return m_coords[i]; }
        inline value_type const &coords(int i) const { return m_coords[i]; }
        inline value_type &x()                          { return m_coords[0]; }
        inline value_type const &x() const              { return m_coords[0]; }
        inline value_type &y()                          { return m_coords[1]; }
        inline value_type const &y() const              { return m_coords[1]; }
        inline value_type &z()                          { return m_coords[2]; }
        inline value_type const &z() const              { return m_coords[2]; }
        inline value_type *extent()                     { return m_extent; }
        inline value_type const *extent() const         { return m_extent; }
        inline value_type &extent(int i)             { return m_extent[i]; }
        inline value_type const &extent(int i) const { return m_extent[i]; }
        inline value_type &u()                          { return m_extent[0]; }
        inline value_type const &u() const              { return m_extent[0]; }
        inline value_type &v()                          { return m_extent[1]; }
        inline value_type const &v() const              { return m_extent[1]; }
        inline value_type &w()                          { return m_extent[2]; }
        inline value_type const &w() const              { return m_extent[2]; }
        inline std::list<particle_type*> &particles() { return m_particles; }
        inline std::list<particle_type*> const &particles() const { return m_particles; }
        inline std::list<box_type*> &children()                        { return m_children; }
        inline std::list<box_type*> const &children() const            { return m_children; }
        inline std::list<box_type*> &colleagues()                      { return m_colleagues; }
        inline std::list<box_type*> const &colleagues() const          { return m_colleagues; }
        inline std::list<box_type*> &lists(int i)                   { return m_lists[i]; }
        inline std::list<box_type*> const &lists(int i) const       { return m_lists[i]; }

        /// @brief Splits the current box in to its children, empty children are ignored.
        ///         If the child box do not contains particles, then it's not created.
        ///
        /// @param particles
        /// @param boxes Contains the physical location for all the boxes.
        /// @return void
        template<typename box_storage_type>
        inline void create_children(box_storage_type &boxes)
        {
            // First, create tentative boxes
            Box box[m_num_children];
            int bits[m_dim];
            for (int i = 0; i < m_num_children; ++i)
            {
                getbits(i,m_dim,bits);
                box[i].u() = 0.5*u();
                box[i].v() = 0.5*v();
                box[i].w() = 0.5*w();
                box[i].x() = x() + bits[0]*box[i].u();
                box[i].y() = y() + bits[1]*box[i].v();
                box[i].z() = z() + bits[2]*box[i].w();
            }
            // Check if boxes are worthy
            for (const_particle_iterator p = m_particles.begin(), end = m_particles.end(); p!= end; ++p)
            {
                box[inchild_id((*p)->position)].particles().push_back(*p);
            }

            // If they are, add them to the list of children
            for (int k = 0; k < m_num_children; ++k)
                if (box[k].particles().size() > 0)
                {
                    box[k].id() = k;
                    box[k].idx() = boxes.size();
                    box[k].level() = m_level + 1;
                    box[k].set_parent(this);
#pragma omp critical
                    {
                        boxes.push_back(box[k]);
                        m_children.push_back(&boxes.back());
                    }
                }

        }


        /**
         * @brief Build list1 of a childless box.  The set consisting of all childless boxes adjacent to this one.  This list empty if is children().size() > 0.
         *
         * @return void
         **/
        inline void buildList1()
        {

            if (m_children.size() == 0)
            {
                if (m_parent)
                {
                    for (box_iterator b = m_parent->colleagues().begin(), end = m_parent->colleagues().end(); b != end; ++b)
                    {
                        box_type *box = *b;
                        if (share_point(box))
                            buildList1(box);
                    }
                    for (box_iterator b = m_parent->children().begin(), end = m_parent->children().end(); b != end; ++b)
                    {
                        box_type *box = *b;
                        if (box != this)
                            buildList1(box);
                    }
                }
            }
            else
            {
                for (box_iterator b = m_children.begin(), end = m_children.end(); b != end; ++b)
                {
#pragma omp task firstprivate(b)
                    {
                        (*b)->buildList1();
                    }
                }
            }
#pragma omp taskwait
        }

        /**
         * @brief Add box to list1 if childless or add childs of box if childless, recursively
         *
         * @param box box to add to list1
         * @return void
         **/
        inline void buildList1(box_type *box)
        {
            if (box->children().size() > 0)
            {
                for (box_iterator b = box->children().begin(), end = box->children().end(); b != end; ++b)
                    if (share_point(*b))
                        buildList1(*b);
            }
            else
            {
                m_lists[0].push_back(box);
            }
        }

        /**
         * @brief Build list2 of a box.  The set consisting of all children of colleagues of m_parent that are well separated from this box.
         *
         * @return void
         **/
        inline void buildList2()
        {
            if (m_parent)
            {
                for (box_iterator b = m_parent->colleagues().begin(), end = m_parent->colleagues().end(); b != end; ++b)
                {
                    box_type *box = *b;
                    for (box_iterator bb = box->children().begin(); bb != box->children().end(); ++bb)
                        if (!share_point(*bb))
                        {
                            m_lists[1].push_back(*bb);
                        }
                }
            }
            for (box_iterator b = m_children.begin(), end = m_children.end(); b != end; ++b)
            {
#pragma omp task firstprivate(b)
                {
                    (*b)->buildList2();
                }
            }

#pragma omp taskwait
        }

        /**
         * @brief Adds box to list3 if satifies criteria
         *
         * @param box
         * @return void
         **/
        inline void buildList3(box_type *box)
        {
            for (box_iterator b = box->children().begin(), end = box->children().end(); b != end; ++b)
                if (share_point(*b))
                    buildList3(*b);
                else
                {
                    m_lists[2].push_back(*b);
                }
        }

        /**
         * @brief Build list3 of a box.  List3 is the set consisting of all descendants of this box's
         *        colleagues that are not adjacent to this one but its parent are adjacent to this one.
         *
         * @return void
         **/
        inline void buildList3()
        {
            if (m_children.size() == 0)
            {
                for (box_iterator b = m_colleagues.begin(), end = m_colleagues.end(); b != end; ++b)
                    buildList3(*b);
            }
            else
            {
                for (box_iterator b = m_children.begin(), end = m_children.end(); b != end; ++b)
                {
#pragma omp task firstprivate(b)
                    {
                        (*b)->buildList3();
                    }
                }
            }
#pragma omp taskwait
        }

        /**
         * @brief Builds list4 of a box.  List4 consist of boxes such that the current box is in its list3.  It assumes that buildList3() has been called.
         *
         * @return void
         **/
        inline void buildList4()
        {
            if (m_children.size() == 0)
            {
                for (box_iterator b = m_lists[2].begin(), end = m_lists[2].end(); b != end; ++b)
                {
                    (*b)->lists(3).push_back(this);
                }
            }
            else
            {
                for (box_iterator b = m_children.begin(), end = m_children.end(); b != end; ++b)
                {
#pragma omp task firstprivate(b)
                    {
                        (*b)->buildList4();
                    }
                }
            }
#pragma omp taskwait
        }

        /**
         * @brief ...
         *
         * @return void
         **/
        inline void buildColleagues()
        {
            // If this box's parent has colleagues, then compare this box with their children only,
            // otherwise, make all children of this box colleague of each other.
            if (m_parent)
            {
                for (box_iterator b = m_parent->colleagues().begin(), end = m_parent->colleagues().end(); b != end; ++b)
                {
                    box_type *box = *b;
                    for (box_iterator bb = box->children().begin(); bb != box->children().end(); ++bb)
                        if (share_point(*bb))
                        {
                            m_colleagues.push_back(*bb);
                        }
                }

            }

            for (box_iterator b = m_children.begin(), end = m_children.end(); b != end;)
            {
                box_type *box = *b;
                ++b;
                for (box_iterator bb = b; bb != end; ++bb)
                {
                    {
                        box->colleagues().push_back(*bb);
                        (*bb)->colleagues().push_back(box);
                    }
                }
            }

            for (box_iterator b = m_children.begin(), end = m_children.end(); b != end; ++b)
            {
#pragma omp task firstprivate(b)
                {
                    (*b)->buildColleagues();
                }
            }
#pragma omp taskwait
        }

        /**
         * @brief Check if the box b share a boundary or point with current box
         *
         * @param b box to query
         * @return bool
         **/
        inline bool share_point(box_type *box)
        {
            bool result = false;
            value_type point[m_dim];
            for (int i = 0; i < m_dim; ++i)
                point[i] = box->coords(i)+box->extent(i);

            result = result || intersect(point,*this);

            for (int i = 0; i < m_dim; ++i)
                point[i] = box->coords(i)-box->extent(i);

            result = result || intersect(point,*this);
            return result;
        }

        /**
         * @brief Check if the point at position is at or within the boundaries of the box
         *
         * @param position Position of the particle
         * @param box
         * @return bool True if particle position is in box, false otherwise
         **/
        bool intersect(const value_type position[m_dim], const box_type &box)
        {
            bool result = true;
            for (int i = 0; i < m_dim; ++i)
                result = result && (box.coords(i)-box.extent(i) <= position[i]) && (position[i] <= box.coords(i) +box.extent(i));
            return result;
        }

        /**
         * @brief Check to what child the point at position belongs to
         *
         * @param position Position of the particle
         * @return box id
         **/
        inline int inchild_id(const value_type position[m_dim])
        {
            int tuple[m_dim];
            for (int k = 0; k < m_dim; ++k)
                tuple[k] = (position[k]-m_coords[k]) < 0 ? 0 : 1;

            return bits2int(m_dim,tuple);
        }

        template<typename output_type>
        friend output_type & operator<< (output_type & o,box_type const & box)
        {
            o.precision(16);
            o << "box idx = " << box.idx() << std::endl;
            o << "   num particles = " << box.m_particles.size() << std::endl;
            o << "   center = " << "[" << box.x() << "," << box.y() << "," << box.z() << "]" << std::endl;
            o << "   extent = " << "[" << box.u() << "," << box.v() << "," << box.w() << "]" << std::endl;
            o << "   label = " << box.level();
            o << " " << box.id();
            if (box.parent()) o << " " << box.parent()->id();
            o << std::endl;

            o << "   childrens = [";
            for (const_box_iterator b = box.children().begin(), end = box.children().end(); b != end; ++b)
                o << (*b)->idx() << " ";
            o << "]" << std::endl;

            o << "   colleagues = [";
            for (const_box_iterator b = box.colleagues().begin(), end = box.colleagues().end(); b != end; ++b)
                o << (*b)->idx() << " ";
            o << "]" << std::endl;

            o << "   list[0] = [";
            for (const_box_iterator b = box.lists(0).begin(), end = box.lists(0).end(); b != end; ++b)
                o << (*b)->idx() << " ";
            o << "]" << std::endl;

            o << "   lists[1] = [";
            for (const_box_iterator b = box.lists(1).begin(), end = box.lists(1).end(); b != end; ++b)
                o << (*b)->idx() << " ";
            o << "]" << std::endl;

            o << "   lists[2] = [";
            for (const_box_iterator b = box.lists(2).begin(), end = box.lists(2).end(); b != end; ++b)
                o << (*b)->idx() << " ";
            o << "]" << std::endl;

            o << "   lists[3] = [";
            for (const_box_iterator b = box.lists(3).begin(), end = box.lists(3).end(); b != end; ++b)
                o << (*b)->idx() << " ";
            o << "]" << std::endl;

            o << "   particles = [";
            for (const_particle_iterator p = box.particles().begin(), end = box.particles().end(); p != end; ++p)
                o << "[" << (*p)->position[0] << "," << (*p)->position[1] << "," << (*p)->position[2] << "];";
            o << "]" << std::endl;

            return o;
        }


};




#endif





