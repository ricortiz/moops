#ifndef ELASTIC_BOUNDARY_HPP
#define ELASTIC_BOUNDARY_HPP
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

#include "particle_system/elastic_system/spring_system.hpp"

/** \class ElasticBoundary
 *  \ingroup ParticleSystem_Module
 *  \brief Handle the spring creation and force computation
 *  \tparam Derived is the derived type, ie the application type or an expression.
 */
template<typename Derived>
class ElasticBoundary : public SpringSystem<Derived>
{
    public:
        typedef SpringSystem<Derived>    	       		base_type;
        typedef typename Traits<Derived>::value_type            value_type;
        typedef typename Traits<Derived>::particle_type         particle_type;
        typedef typename base_type::spring_type       		spring_type;
        typedef typename base_type::spring_iterator   		spring_iterator;
        typedef typename base_type::spring_lut_type             spring_map_type;

    public:

        inline Derived &derived()
        {
            return *static_cast<Derived*>(this);
        }
        
        inline void setSprings(std::vector<size_t> &col_ptr, std::vector<size_t> &col_idx, std::vector<value_type> &strength)
        {
            particle_type *particles = derived().particles();
            for(size_t p = 0; p < col_ptr.size() - 1; ++p)
                for(size_t i = col_ptr[p], end = col_ptr[p + 1]; i < end; ++i)
                    if(!this->existSpring(&particles[p], &particles[col_idx[i]]))
                    {
                        spring_iterator s = this->addSpring(&particles[p], &particles[col_idx[i]], strength[i]);
                        s->getAidx() = 3 * p;
                        s->getBidx() = 3 * col_idx[i];
                    }
        }
                
        inline void computeForces()
        {
            for (spring_iterator s = this->springs_begin(), end = this->springs_end(); s != end; ++s)
                s->apply();
        }

        inline void computeForces(const value_type *x)
        {
            for (spring_iterator s = this->springs_begin(), end = this->springs_end(); s != end; ++s)
            {
                size_t i = s->getAidx();
                size_t j = s->getBidx();
                const value_type *x1 = &x[i];
                const value_type *x2 = &x[j];
                s->apply(x1, x2);
            }
        }

        inline void computeForces(const value_type *x, value_type *f)
        {
            for (spring_iterator s = this->springs_begin(), end = this->springs_end(); s != end; ++s)
            {
                size_t i = s->getAidx();
                size_t j = s->getBidx();
                const value_type *x1 = &x[i];
                const value_type *x2 = &x[j];
                value_type *f1 = &f[i];
                value_type *f2 = &f[j];
                s->apply(x1, x2, f1, f2);
            }
        }

};


#endif
