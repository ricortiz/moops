#ifndef GMRES_RAW_HPP
#define GMRES_RAW_HPP

#include<algorithm>

#include "krylov_storage.hpp"
#include "krylov_base.hpp"

template<typename value_type, size_t system_size, size_t krylov_space_max_dim>
class GeneralizedMinimalResidualMethod : public KrylovBase<GeneralizedMinimalResidualMethod<value_type, system_size,krylov_space_max_dim> >
{
    protected:
        krylov_storage<value_type,system_size,krylov_space_max_dim> m_storage;

    private:
        enum
        {
            k_max = krylov_space_max_dim
        };

    public:

        inline value_type &H ( size_t i,size_t j ) { return m_storage.H ( i,j ); }
        inline value_type *H () { return m_storage.H (); }
        inline value_type *residual() { return m_storage.residual(); }
        inline value_type &residual ( size_t i ) { return m_storage.residual ( i ); }
        inline value_type &s ( size_t i ) { return m_storage.s ( i ); }
        inline value_type &c ( size_t i ) { return m_storage.c ( i ); }
        inline value_type &g ( size_t i ) { return m_storage.g ( i ); }
        inline value_type *g ( ) { return m_storage.g(); }
        inline value_type *v ( size_t i ) { return m_storage.v ( i ); }

        /**
        * \brief Gneralized Minimal Residual Method.  Solves the linear system:  F'(xc)*x = -f0.
        *
        * \param F This is a linear operator that acts on x.
        * \param b Right hand side vector of the linear system.
        * \param x Initial guess and solution vector
        * \param stats This is an optional vector to collect statistics of the method.  Defaults to 0.
        * \return norm of the residual
        **/
        template<typename operator_type>
        inline value_type operator() ( operator_type &F, const value_type *b, value_type *x, value_type errtol = 1e-6, value_type *stats = 0 )
        {
            F ( x,residual() );
            std::transform ( b,b+system_size,residual(),residual(),std::plus<value_type>() );
            value_type rho = norm ( residual() );
            errtol *= rho;

            /// Test for early termination
            if ( rho < errtol )
                return rho;

            std::fill ( g(),g() +k_max+1,value_type ( 0 ) );
            g ( 0 ) = rho;

            for ( size_t i = 0 ; i < system_size; ++i )
                v ( 0 ) [i] = -residual ( i ) / rho;
            /// Start GMRES iteration
            unsigned int k;
            for ( k = 0; k < k_max && rho > errtol; ++k )
            {
                /// Evaluate linear operator
                F ( v ( k ),v ( k+1 ) );

                /// Apply Arnoldi's method with modified Gram-Schmidt orthogonalization to the vectors v
                value_type ld = this->arnoldi ( k );

                /// Perform Givens rotations in order to eliminate the H(k+1,k) entry from the
                /// Hessenberg matrix.
                for ( size_t i = 0; i < k; ++i )
                    apply_rotation ( H ( i,k ), H ( i+1,k ), c ( i ), s ( i ) );
                get_rotation ( H ( k,k ), ld, c ( k ), s ( k ) );
                apply_rotation ( H ( k,k ), ld, c ( k ), s ( k ) );
                apply_rotation ( g ( k ), g ( k+1 ), c ( k ), s ( k ) );

                assert ( ld < 1e-16 || !"rotation failed: Lower diagonal is non-zero" );
                ///  The norm of the residual is the last entry in g
                rho = std::fabs ( g ( k+1 ) );
            }

            /// Solve the least square problem by solving a upper triangular linear system
            this->back_solve ( k );
            /// Update the solution
            for ( int i = 0; i < k; i++ )
                for ( size_t j = 0; j < system_size; ++j )
                    x[j] += g ( i ) * v ( i ) [j];

            return rho;
        }
};

template<typename _value_type, size_t _system_size, size_t _krylov_space_max_dim>
struct krylov_traits<GeneralizedMinimalResidualMethod<_value_type, _system_size, _krylov_space_max_dim> >
{
    typedef _value_type value_type;
    enum
    {
        krylov_space_max_dim = _krylov_space_max_dim,
        system_size = _system_size
    };
};

#endif
