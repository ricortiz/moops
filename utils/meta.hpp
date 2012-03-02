#ifndef SDC_META_HPP
#define SDC_META_HPP
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
template<bool conditional, typename Then, typename Else>
struct IF
{
    typedef Then type;
};
template<typename Then, typename Else>
struct IF<false, Then, Else>
{
    typedef Else type;
};

template<bool conditional, typename Then, typename Else>
struct IFP : public IF<conditional, Then, Else>::type {};

template<bool condition>
struct StaticAssert;

template<>
struct StaticAssert<true> { enum {WRONG_QUADRATURE_CHOSEN}; };

#define STATIC_ASSERT(CONDITION,MSG) \
    if (StaticAssert<bool(CONDITION)>::MSG) {}

#define STATIC_ASSERT_QUAD(TYPE) \
    STATIC_ASSERT((TYPE::quadrature == gauss_lobatto \
                || TYPE::quadrature == gauss_radau \
                || TYPE::quadrature == clenshaw_curtis), \
                WRONG_TYPE_OF_QUADRATURE)



#endif
