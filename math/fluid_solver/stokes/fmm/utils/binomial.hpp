#ifndef FMM_BINOMIAL_HPP
#define FMM_BINOMIAL_HPP
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
#include<sstream>
#include<cmath>
#include<stdexcept>
#include<typeinfo>

#include"factorial.hpp"

template <class T>
T binomial_coefficient(unsigned n, unsigned k)
{
    static const char* function = "binomial<%1%>(unsigned, unsigned)";
    if (k > n)
        return 0;
    T result;
    if ((k == 0) || (k == n))
        return 1;
    if ((k == 1) || (k == n-1))
        return n;

    if (n <= max_factorial<T>::value)
    {
        // Use fast table lookup:
        result = factorial<T>(n);
        result /= (factorial<T>(n-k)*factorial<T>(k));
    }
    else 
    {
        std::stringstream s;
        s << function << ": The binomial coefficient is undefined for n > " << max_factorial<T>::value;
        throw std::domain_error(s.str());
    }
    return result;
    // convert to nearest integer:
//     return std::ceil(result - 0.5f);
}


#endif
