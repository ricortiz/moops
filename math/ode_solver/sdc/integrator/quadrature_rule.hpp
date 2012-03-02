#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP
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
#include "utils/meta.hpp"
#include "gauss_lobatto.hpp"
#include "gauss_radau.hpp"
#include "clenshaw_curtis.hpp"

enum quadrature_type { gauss_lobatto, gauss_radau, clenshaw_curtis, other};

template<typename value_type, quadrature_type quadrature, int sdc_nodes, int multirate_nodes>
struct QuadratureRule;

template<typename value_type, quadrature_type quadrature, int sdc_nodes>
struct QuadratureRule<value_type, quadrature, sdc_nodes, 0>
{
    typedef IFP < quadrature == gauss_lobatto, GaussLobatto<value_type, sdc_nodes>,
            IFP < quadrature == gauss_radau, GaussRadau<value_type, sdc_nodes>,
            IFP < quadrature == clenshaw_curtis, ClenshawCurtis<value_type, sdc_nodes>, StaticAssert<true> > > > type;
};

#endif
