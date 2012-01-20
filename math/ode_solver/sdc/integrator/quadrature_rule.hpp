#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include "../utils/meta.hpp"
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
