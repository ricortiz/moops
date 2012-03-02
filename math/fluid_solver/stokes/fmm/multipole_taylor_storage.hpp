#ifndef MULTIPOLE_TAYLOR_STORAGE_HPP
#define MULTIPOLE_TAYLOR_STORAGE_HPP
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
#include "math/fluid_solver/stokes/fmm/utils/meta.hpp"
#include "math/fluid_solver/stokes/fmm/utils/factorial.hpp"
#include "math/fluid_solver/stokes/fmm/utils/binomial.hpp"

template<typename T, int precision, int num_coefficients  >
struct multipole_taylor_arrays
{
    T* phi;
    T* psi;
    T ijk_factorial[num_coefficients];
    T ijk_binomial[num_coefficients];

    T binomial[precision+1][precision+1];
    int tetra_two[precision+3];
    int tetra_three[precision+3];

};

/** \internal
 *
 * \class multipole_taylor_storage
 *
 * \brief Stores the data of the multipole_taylor
 *
 * This class stores the data
 *
 */
template<typename T, int precision, int num_coefficients>
class multipole_taylor_storage;


template<typename T, int precision, int _num_coefficients = Binomial<precision+3,3>::value>
class multipole_taylor_storage
{
    private:
        multipole_taylor_arrays<T,precision,_num_coefficients> m_data;
        size_t m_coeff_size;

    public:
        explicit multipole_taylor_storage(size_t num_boxes) : m_coeff_size(num_boxes*_num_coefficients)
        {
            m_data.phi = new T[4*m_coeff_size];
            m_data.psi = new T[4*m_coeff_size];

            for (size_t i = 0, end = 4*m_coeff_size; i < end; ++i)
            {
                m_data.phi[i] = 0;
                m_data.psi[i] = 0;
            }

            for (int i = 0, idx = 0; i <= precision; ++i)
                for (int j = 0; j <= precision-i; ++j)
                {
                    for (int k = 0; k <= precision-i-j; ++k, ++idx)
                    {
                        ijk_factorial(idx) = T(1) / (factorial<T> (i) *factorial<T> (j) *factorial<T> (k));
                        ijk_binomial(idx) = factorial<T> (i+j+k) *ijk_factorial(idx);
                    }
                    binomial(i,j) = binomial_coefficient<T> (i,j);
                }

            for (int i = 0; i <= precision+2 ; ++i)
            {
                tetra_three(i) = num_coefficients-binomial_coefficient<int> (precision-i+2,3);
                tetra_two(i) = binomial_coefficient<int> (precision-i+2,2);
            }
        }

        ~multipole_taylor_storage()
        {
            delete [] m_data.phi;
            delete [] m_data.psi;
        }
        inline void swap(multipole_taylor_storage& other) { std::swap(m_data,other.m_data); }

        inline const T *phi(int i, int j) const { return m_data.phi+i*num_coefficients+j*m_coeff_size; }
        inline T *phi(int i, int j)             { return m_data.phi+i*num_coefficients+j*m_coeff_size; }

        inline const T *psi(int i, int j) const { return m_data.psi+i*num_coefficients+j*m_coeff_size; }
        inline T *psi(int i, int j)             { return m_data.psi+i*num_coefficients+j*m_coeff_size; }

        inline const T &ijk_factorial(int i) const { return m_data.ijk_factorial[i]; }
        inline T &ijk_factorial(int i)             { return m_data.ijk_factorial[i]; }
        
        inline const T *ijk_factorial() const { return m_data.ijk_factorial; }
        inline T *ijk_factorial()             { return m_data.ijk_factorial; }

        inline const T &ijk_binomial(int i) const { return m_data.ijk_binomial[i]; }
        inline T &ijk_binomial(int i)             { return m_data.ijk_binomial[i]; }

        inline const T &binomial(int i,int j) const { return m_data.binomial[i][j]; }
        inline T &binomial(int i,int j)             { return m_data.binomial[i][j]; }

        inline const int &tetra_two(int i) const { return m_data.tetra_two[i]; }
        inline int &tetra_two(int i)             { return m_data.tetra_two[i]; }

        inline const int &tetra_three(int i) const { return m_data.tetra_three[i]; }
        inline int &tetra_three(int i)             { return m_data.tetra_three[i]; }
    
        inline size_t coeff_size() { return 4*m_coeff_size; }

        enum
        {
            num_coefficients = _num_coefficients
        };
};


#endif
