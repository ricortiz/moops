#ifndef NEWTON_STORAGE_HPP
#define NEWTON_STORAGE_HPP
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
/** \internal
 *
 */
template<typename T>
struct newton_arrays
{
    T *dx;
    T *f;
};

/** \internal
 *
 *  \class newton_storage
 *
 * \brief Stores the data used in the newton method
 *
 *
 */
template<typename T>
class newton_storage;

template<typename T>
class newton_storage
{
        newton_arrays<T> m_data;
    public:
        inline explicit newton_storage (size_t size)
        {
            m_data.dx = new T[size];
            m_data.f = new T[size];
            clear(size);
        }
        void clear(size_t size)
        {
            std::fill(m_data.dx,m_data.dx+size,0.0);
            std::fill(m_data.f,m_data.f+size,0.0);
        }
        ~newton_storage()
        {
            delete [] m_data.dx;
            delete [] m_data.f;
        }
        inline void swap ( newton_storage& other ) { std::swap ( m_data,other.m_data ); }

        inline T *dx() { return m_data.dx; }
        inline T &dx(size_t i) { return * ( m_data.dx+i); }
        inline T *f ( ) { return m_data.f; }
};


#endif
