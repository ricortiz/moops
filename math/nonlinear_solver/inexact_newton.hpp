#ifndef INEXACT_NEWTON_RAW_HPP
#define INEXACT_NEWTON_RAW_HPP

#include<cmath>
#include<limits>
#include<vector>

#include "math/nonlinear_solver/newton_storage.hpp"
#include "math/nonlinear_solver/newton_base.hpp"
#include "math/linear_solver/krylov/generalized_minimal_residual_method_raw.hpp"


/**
 * \brief Matrix free Newton-Krylov gloabaly convergent nonlinear solver.
 * Inexact-Newton-Armijo iteration.
 * Eisenstat-Walker forcing term.
 * Parabolic line search via three point interpolation.
 * Computes Jacobian using finite differences.
 **/
template< typename value_type, int k_max, int k_restart, typename krylov_method_type = GeneralizedMinimalResidualMethod<value_type,k_max> >
class InexactNewtonMethod;

template< typename value_type,int k_max,int k_restart >
class InexactNewtonMethod<value_type,k_max,k_restart,GeneralizedMinimalResidualMethod<value_type,k_max> >
    : public NewtonBase<InexactNewtonMethod<value_type,k_max,k_restart,GeneralizedMinimalResidualMethod<value_type,k_max> > >
{
    private:
        size_t m_system_size;
        newton_storage<value_type> 	m_storage;
        GeneralizedMinimalResidualMethod<value_type,k_max> m_gmres;
        value_type 			m_gamma;
        value_type 			m_etamax;
	
    public:

        InexactNewtonMethod(size_t system_size) :
                m_system_size(system_size),
                m_storage(system_size),
                m_gmres(system_size),
                m_gamma ( value_type ( .9 ) ),
                m_etamax ( value_type ( .9 ) )
        {};

        inline value_type *dx() { return m_storage.dx(); }
        inline value_type *f() { return m_storage.f(); }
        inline value_type *ft() { return m_storage.ft(); }
        inline value_type *xt() { return m_storage.xt(); }
        inline value_type &xt ( size_t i ) { return m_storage.xt ()[i]; }
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
        void operator() ( operator_type &F, value_type *x, value_type atol, value_type rtol, std::vector<value_type> *stats = 0 )
        {
            typedef directional_derivative<operator_type,value_type>   jacobian_operator_type;
            jacobian_operator_type jacobian ( F,x,f(),m_system_size);
            /// Evaluate F at initial iterate and compute the stop tolerance

            F ( x,f() );
            value_type fnrm = norm ( f() );
            value_type fnrmo = value_type ( 1 );
            value_type stop_tol = atol + rtol * fnrm;
            value_type gmres_tol = m_etamax;

            unsigned int itc = 0;

            /// Main iteration loop
            while ( fnrm > stop_tol )
            {
                /// Compute ratio of succesive residual norms and iteration counter
                itc++;
                value_type rat = fnrm / fnrmo;
                fnrmo = fnrm;
                std::fill ( dx(),dx() +m_system_size,value_type ( 0 ) );

                value_type k_err = std::numeric_limits<value_type>::infinity();

                /// Solve for descend direction using GMRES
                unsigned int k_it = 0;
                while ( k_err > gmres_tol*fnrm && k_it++ < k_restart && fnrm != 0 )
                    k_err = m_gmres ( jacobian, f(), dx(), gmres_tol );

                /// Start Armijo line search
                value_type lambda[3]    = {1., 1., 1.};
                value_type fnorm[2]     = {0., 0.};
                value_type fnorm_sqr[3] = {0., 0., 0.};

                for ( size_t i = 0; i < m_system_size; ++i )
                    xt ( i ) = x[i] + lambda[0] * dx ( i );

                F ( xt(), ft() );
                fnorm[0] = norm ( f() );
                fnorm[1] = norm ( ft() );
                fnorm_sqr[0] = fnorm[0] * fnorm[0];
                fnorm_sqr[1] = fnorm[1] * fnorm[1];
                fnorm_sqr[2] = fnorm_sqr[1];
                armijo(F,x,lambda,fnorm,fnorm_sqr);

                std::copy ( xt(),xt()+m_system_size,x );
                std::copy ( ft(),ft()+m_system_size,f() );
                /// End of Armijo line search.
                fnrmo = fnrm;
                fnrm = norm ( f() );
                rat = fnrm / fnrmo;

                /// Adjust eta as per Eisenstat-Walker.
                if ( m_etamax > 0 )
                {
                    value_type etaold = m_etamax;
                    value_type etanew = m_gamma * rat * rat;
                    if ( m_gamma*etaold*etaold > .1 )
                        etanew = std::max ( etanew, m_gamma * etaold * etaold );
                    gmres_tol = std::max ( std::min ( etanew, etaold ), value_type ( 0.5 ) * stop_tol / fnrm );
                }
            }
            if ( fnrm > stop_tol )
                std::cout << "Newton method failed to converge to desired accuracy." << std::endl;
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

