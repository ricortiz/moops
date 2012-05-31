#ifndef INEXACT_NEWTON_HPP
#define INEXACT_NEWTON_HPP
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
#include<limits>
#include<vector>
#include<map>
#include<sstream>

#include "math/nonlinear_solver/newton_storage.hpp"
#include "math/nonlinear_solver/newton_base.hpp"
#include "math/linear_solver/krylov/generalized_minimal_residual_method.hpp"


/**
 * \brief Matrix free Newton-Krylov gloabaly convergent nonlinear solver.
 * Inexact-Newton-Armijo iteration.
 * Eisenstat-Walker forcing term.
 * Parabolic line search via three point interpolation.
 * Computes Jacobian using finite differences.
 **/
template< typename value_type,int k_max,int k_restart, typename linear_solver_type = GeneralizedMinimalResidualMethod<value_type,k_max> >
class InexactNewtonMethod: public NewtonBase<InexactNewtonMethod<value_type,k_max,k_restart,linear_solver_type> >
{
    private:
        size_t m_system_size;
        size_t m_maxitc;
        newton_storage<value_type> 	m_storage;
        linear_solver_type m_gmres;
        value_type 			m_gamma;
        value_type 			m_forcing_term;
        
    public:

        InexactNewtonMethod(size_t system_size) :
                m_system_size(system_size),
                m_maxitc(10),
                m_storage(system_size),
                m_gmres(system_size),
                m_gamma ( value_type ( .9 ) ),
                m_forcing_term ( value_type ( .9 ) )
        {};

        inline value_type *dx() { return m_storage.dx(); }
        inline value_type *f() { return m_storage.f(); }
        inline value_type &dx ( size_t i ) { return m_storage.dx ()[i]; }
        inline size_t system_size() { return m_system_size; }

        /**
        * \brief Newton iterations routine.  It solves F(x) = 0.
        *
        * \param F Non-Linear Function.  This is actually a functor (function object).
        * \param x Initial guess and result vector.
        * \param atol Absolute tolerance.
        * \param rtol Relative tolerance.
        * \param stats Vector to collect statistics about the method. Defaults to 0.
        *
        **/
        template< typename operator_type>
        int operator() ( operator_type &F, value_type *x, value_type atol = 1e-3, value_type rtol = 1e-4, std::map<std::string,std::vector<value_type> > *stats = 0 )
        {
            typedef directional_derivative<operator_type,value_type>   jacobian_operator_type;
            jacobian_operator_type jacobian ( F,x,f(),m_system_size);
            /// Evaluate F at initial iterate and compute the stop tolerance

            F ( x,f() );
            value_type fnrm = this->norm ( f(),m_system_size );
            value_type fnrmo = value_type ( 1 );
            value_type stop_tol = atol + rtol * fnrm;
            value_type gmres_tol = m_forcing_term;
            if(stats)
            {
                std::map<std::string,std::vector<value_type> > &s = *stats;
                std::string iterate_key = "newton_iterate_0";
                std::string norm_key = "newton_norm";
                std::string fn_eval_newton_key = "newton_fn_eval";
                std::string armijo_it_key = "armijo_fn_evals";

                s[norm_key].push_back(fnrm);
                s[fn_eval_newton_key].push_back(0);
                s[armijo_it_key].push_back(0);
                s[iterate_key].resize(m_system_size);
                std::copy(x,x+m_system_size,s[iterate_key].begin());
            }

            unsigned int itc = 0;

            
            while ( fnrm > stop_tol && itc < m_maxitc )                         // Main iteration loop
            {
                itc++;                                                          // Newton iteration counter
                value_type rat = fnrm / fnrmo;                                  // Compute ratio of succesive residual norms and iteration counter
                fnrmo = fnrm;                                                   // Store old function norm
                std::fill ( dx(),dx() +m_system_size,value_type ( 0 ) );        // Set initial Krylov iterate to zero
                value_type k_err = std::numeric_limits<value_type>::infinity(); // Define initial error to infinity                
                unsigned int k_it = 0;
                while ( k_err > gmres_tol*fnrm && k_it++ < k_restart && fnrm != 0 ) // Solve for descend direction using GMRES
                    k_err = m_gmres ( jacobian, f(), dx(), gmres_tol, stats );

                std::transform(x, x + m_system_size,dx(),x,std::minus<value_type>()); // Update x: x = x + -dx.  dx = steppest descent direction
                                
                value_type lambda[3]    = {1., 1., 1.};                               // Start Armijo line search
                value_type fnorm[2]     = {0., 0.};
                value_type fnorm_sqr[3] = {0., 0., 0.};
                
                fnorm[0] = this->norm ( f(), m_system_size );
                F ( x, f() );
                fnorm[1] = this->norm ( f(), m_system_size );
                fnorm_sqr[0] = fnorm[0] * fnorm[0];
                fnorm_sqr[1] = fnorm[1] * fnorm[1];
                fnorm_sqr[2] = fnorm_sqr[1];

                if ( !this->armijo(F,x,lambda,fnorm,fnorm_sqr,stats) )
                {
                    if (stats)
                    {
                        std::map<std::string,std::vector<value_type> > &s = *stats;
                        std::stringstream iterate_key;
                        iterate_key << "newton_iterate_" << itc;
                        s[iterate_key.str()].resize(m_system_size);
                        std::copy(x,x+m_system_size,s[iterate_key.str()].begin());
                    }
                    return 2;
                }
                ///< End of Armijo line search.
                
                fnrmo = fnrm;
                fnrm = fnorm[1];
                rat = fnrm / fnrmo;
                if(stats)
                {
                    std::map<std::string,std::vector<value_type> > &s = *stats;
                    std::stringstream iterate_key;
                    iterate_key << "newton_iterate_" << itc;
                    std::string norm_key = "newton_norm";
                    std::string fn_eval_newton_key = "newton_fn_eval";
                    std::string armijo_it_key = "armijo_fn_evals";
                    std::string fn_eval_gmres_key = "gmres_fn_eval";
                    s[iterate_key.str()].resize(m_system_size);
                    std::copy(x,x+m_system_size,s[iterate_key.str()].begin());
                    s[norm_key].push_back(fnrm);
                    value_type total_fn_evals = s[fn_eval_newton_key].back() // Function evaluations from previous iterations
                            + s[fn_eval_gmres_key].back() // Function evaluations from GMRES
                            + s[armijo_it_key].back() + 1; // Function evaluations from the armijo iterations
                    if(itc == 1) total_fn_evals++;
                    s[fn_eval_newton_key].push_back(total_fn_evals);
                }
                /// Adjust eta as per Eisenstat-Walker.
                if ( m_forcing_term > 0 )
                {
                    value_type etaold = m_forcing_term;
                    value_type etanew = m_gamma * rat * rat;
                    if ( m_gamma*etaold*etaold > .1 )
                        etanew = std::max ( etanew, m_gamma * etaold * etaold );
                    gmres_tol = std::max ( std::min ( etanew, etaold ), value_type ( 0.5 ) * stop_tol / fnrm );
                }
            }
            
            if ( fnrm > stop_tol )
            {
                std::cout << "Newton method failed to converge to desired accuracy. Function Norm = " << fnrm << std::endl;
                return 1;
            }
            return 0;
        }

};

template<typename _value_type,int _k_max, int _k_restart, typename _krylov_method_type>
struct newton_traits<InexactNewtonMethod<_value_type,_k_max,_k_restart,_krylov_method_type> >
{
    typedef _value_type value_type;
    enum
    {
        k_max = _k_max,
        k_restart = _k_restart
    };
};




#endif

