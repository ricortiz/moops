#ifndef PLUMMER_DISTRIBUTION_HPP
#define PLUMMER_DISTRIBUTION_HPP
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
#include<cstdlib>

template<typename value_type>
void genereatePoint(value_type rad,value_type *p)
{
    p[0] = 0;
    p[1] = 0;
    p[2] = 0;
    bool tooBig = true;
    value_type rsq = 0;
    while (rsq <= 1)
    {
        rsq = 0;
        p[0] = 2*rand()/((value_type)RAND_MAX+1) - 1;
        rsq = rsq + p[0] * p[0];
        p[1] = 2*rand()/((value_type)RAND_MAX+1) - 1;
        rsq = rsq + p[1] * p[1];
        p[2] = 2*rand()/((value_type)RAND_MAX+1) - 1;
        rsq = rsq + p[2] * p[2];
    }
    
    p[0] = rad * p[0]/std::sqrt(rsq);
    p[1] = rad * p[1]/std::sqrt(rsq);
    p[2] = rad * p[2]/std::sqrt(rsq);
    
}

template<typename value_type>
value_type genereatePoint(value_type rad)
{
    value_type p = 0;
    bool tooBig = true;
    value_type rsq = 0;
    while (rsq <= 1)
    {
        rsq = 0;
        p = 2*rand()/((value_type)RAND_MAX+1) - 1;
        rsq = rsq + p * p;
    }
    
    return p;
    
}


template<typename value_type>
void plummerDistribution(value_type position[], size_t num_particles)
{
    //used in distribution creation
    value_type rsc = 3 * M_PI/16;
    
    value_type h = 1.0/(num_particles-1);
    for (int i = 0, idx = 0; i < num_particles; i++, idx+=3)
    {
        value_type temp = 1.0/std::pow(0.999*rand()/((value_type)RAND_MAX+1),1.5) - 1.0;
        value_type ri = 1.0/sqrt(temp);
        
        genereatePoint(rsc*ri,&position[idx]);
    }
    
}

#endif
