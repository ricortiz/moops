#ifndef OCTREE_STORAGE_HPP
#define OCTREE_STORAGE_HPP
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
template<typename box_type, typename array_type>
struct octree_arrays
{
    array_type       boxes;
};


/** \internal
 *
 * \class octree_storage
 *
 * \brief Stores the data of the octree
 *
 * This class stores the data
 *
 */
template<typename box_type, typename array_type>
class octree_storage;


template<typename box_type, typename array_type>
class octree_storage
{
    private:
        octree_arrays<box_type,array_type> m_data;

    public:
        explicit octree_storage()
        {
            // This creates the root of the tree
            m_data.boxes.push_back ( box_type() );
        }
        ~octree_storage() {}
        inline void swap ( octree_storage& other ) { std::swap ( m_data,other.m_data ); }

        inline const box_type &root() const { return m_data.boxes.front(); }
        inline box_type &root()             { return m_data.boxes.front(); }

        inline const array_type &boxes() const { return m_data.boxes; }
        inline array_type       &boxes() { return m_data.boxes; }
};

#endif


