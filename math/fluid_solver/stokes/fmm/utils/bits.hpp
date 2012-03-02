#ifndef BITS_HPP
#define BITS_HPP
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

/// @brief Extracts the n least significant bits of x and puts them in
///         the array "bits" replacing 0 by -1
/// bits[0] = most significant bit, bits[1] = next bit, etc.
void getbits ( int x, int n, int bits[] )
{
    int pattern;

    pattern = 1 << (n-1);
    for ( int i = 0; i < n; ++i, pattern >>= 1 )
        bits[i] = ( x & pattern ) ? 1 : -1;

}


/// @brief Convert array of bits to an integer
/// bits[] has same format as in getbits()
/// (i.e., bits[0] is m.s.b. and bits[n-1] is l.s.b.)
int bits2int ( int n, int bits[] )
{
    int x;
    int pattern = 1 << (n-1);
    x = 0;
    for ( int i = 0; i < n; i++, pattern >>= 1 )
        if ( bits[i] == 1 )
            x |= pattern;

    return x;
}

#endif
