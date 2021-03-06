#ifndef NEWTON_BASE_HPP
#define NEWTON_BASE_HPP
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
#include<iterator>
#include<cassert>

template<typename T>
class newton_traits;

template<typename Derived>
class NewtonBase
{
    public:
        typedef typename newton_traits<Derived>::value_type value_type;

    private:
        value_type m_alpha;
        value_type m_sigma0; //< safeguarding bounds for the linesearch
        value_type m_sigma1; //< safeguarding bounds for the linesearch
        size_t m_maxarm;

    public:
        NewtonBase() : m_alpha ( value_type ( 1e-4 ) ), m_sigma0 ( value_type ( .1 ) ), m_sigma1 ( value_type ( .5 ) ), m_maxarm ( 20 )
        {}

        /**
         * @brief This function returns an instance of the derived type.
         *
         * @return Derived&
         **/
        Derived &derived()
        {
            return *static_cast<Derived*> ( this );
        }

        template< typename operator_type, typename vector_type>
        inline bool armijo ( operator_type &F, value_type *x, value_type lambda[3], value_type fnorm[2], value_type fnorm_sqr[3], std::map<std::string,vector_type> *stats = 0 )
        {
            size_t system_size = derived().system_size();
            unsigned int iarm = 0; // Armijo iteration counter
            while ( fnorm[1] >= ( value_type ( 1 ) - m_alpha * lambda[0] ) *fnorm[0] )
            {
                /// Apply three point parabolic model
                if ( iarm == 0 )
                    lambda[0] *= m_sigma1;
                else
                    lambda[0] = parab3p ( lambda[1], lambda[2], fnorm_sqr[0], fnorm_sqr[1], fnorm_sqr[2] );

                /// Update x and keep books on lambda
                for ( size_t i = 0; i < system_size; ++i )
                    x[i] -= lambda[0] * derived().dx ( i );
                lambda[2] = lambda[1];
                lambda[1] = lambda[0];
                /// Keeps books on function norms
                F ( x, derived().f() );

                fnorm[1] = norm ( derived().f(), system_size );
                fnorm_sqr[2] = fnorm_sqr[1];
                fnorm_sqr[1] = fnorm[1] * fnorm[1];
                iarm++;

                if ( iarm > m_maxarm )
                {
                    if(stats)
                    {
                        std::map<std::string,vector_type> &s = *stats;
                        std::string armijo_it_key = "armijo_fn_evals";
                        s[armijo_it_key].push_back(iarm);
                    }
                    return false;
                }
            }
            if(stats)
            {
                std::map<std::string,vector_type> &s = *stats;
                std::string armijo_it_key = "armijo_fn_evals";
                s[armijo_it_key].push_back(iarm);
            }
            return true;
        }

        /**
         * \brief Three-point safeguarded parabolic model for a line search.
         *
         * \param lamc Current steplength
         * \param lamm Previous steplength
         * \param ff0 Value of |F(xc)|^2
         * \param ffc Value of |F(xc + lamc*d)|^2
         * \param ffm Value of |F(xc + lamm*d)|^2
         * \return New value of lambda given by the parabolic model.
         *
         **/
        inline value_type parab3p ( value_type lamc, value_type lamm, value_type ff0, value_type ffc, value_type ffm )
        {
            /// Compute coefficients of interpolation polynomial
            ///  p(lambda) = ff0 + (c1*lambda + c2*lambda^2)/d1
            ///  d1 = (lamc - lamm)*lamc*lamm < 0
            ///    if c2 > 0, then we have a concave up curvature and defaults to
            ///    lambda = m_sigma1*lambda
            value_type c2 = lamm * ( ffc - ff0 ) - lamc * ( ffm - ff0 );
            if ( c2 >= 0 )
                return m_sigma0 * lamc;

            value_type lambda = 0;
            value_type c1 = lamc * lamc * ( ffm - ff0 ) - lamm * lamm * ( ffc - ff0 );
            lambda = -.5 * c1 / c2;

            if ( lambda < m_sigma0 * lamc )
                lambda = m_sigma0 * lamc;
            if ( lambda > m_sigma1 * lamc )
                lambda = m_sigma1 * lamc;

            return lambda;

        }

        inline value_type dot ( const value_type *x, const value_type *y, size_t system_size )
        {
            value_type s = value_type ( 0 );
            for ( size_t i = 0; i < system_size; ++i )
                s += x[i] * y[i];
            assert(!isnan(s));
            return s;
        }
        inline value_type norm ( const value_type *x, size_t system_size )
        {
            return std::sqrt ( dot ( x, x, system_size ) );
        }

};

/**
 * \brief Finite difference directional derivative.
 *   Computes the directional derivative of the nonlinear function F as an approximation to
 *    F'(xc)*x.
 *
 * \param x Point.
 * \param w Direction vector.
 * \param F Non-linear function.
 * \param f0 F(x)
 * \param z Resultant derivative.
 * \param eps Difference increment. Defaults to value_type(1e-7).
 *
 **/
template<typename operator_type, typename value_type>
struct directional_derivative
{

    operator_type &m_F;
    const value_type *m_x;
    const value_type *m_f0;
    size_t m_system_size;

    directional_derivative ( operator_type &F, const value_type *x, const value_type *f0, size_t system_size ) : m_F ( F ), m_x ( x ), m_f0 ( f0 ), m_system_size ( system_size ) {}

    inline void operator() ( const value_type *w, value_type *DF, value_type eps = value_type ( 1e-9 ) )
    {
        value_type wnorm = norm ( w, m_system_size );
        if ( wnorm == 0 )
        {
            std::fill ( DF, DF + m_system_size, value_type ( 0 ) );
            return;
        }

        /// Scale the step
        value_type xs = dot ( m_x, w, m_system_size ) / wnorm;

        if ( xs != value_type ( 0 ) )
            eps = eps * std::max ( std::fabs ( xs ), value_type ( 1 ) ) * ( xs < 0 ? value_type ( -1 ) : value_type ( 1 ) );

        /// Scale the difference increment
        eps /= wnorm;

        std::vector<value_type> f1(m_system_size);
        std::vector<value_type> x1(m_system_size);
        for ( size_t i = 0; i < m_system_size; ++i )
            x1[i] = m_x[i] + eps * w[i];
        m_F ( &x1[0], &f1[0] );
        value_type h = value_type ( 1 ) / eps;
        for ( size_t i = 0; i < m_system_size; ++i )
            DF[i] = ( f1[i] - m_f0[i] ) * h;
    }
    inline value_type dot ( const value_type *x, const value_type *y, size_t system_size )
    {
        value_type s = value_type ( 0 );
        for ( size_t i = 0; i < system_size; ++i )
            s += x[i] * y[i];
        return s;
    }
    inline value_type norm ( const value_type *x, size_t system_size )
    {
        return std::sqrt ( dot ( x, x, system_size ) );
    }
};

#endif
