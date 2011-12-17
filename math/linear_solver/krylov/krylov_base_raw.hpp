#ifndef KRYLOV_BASE_RAW_HPP
#define KRYLOV_BASE_RAW_HPP

#include <cmath>
#include<cassert>

template<typename T> struct krylov_traits;

template<typename Derived>
class KrylovBase
{
    protected:
        typedef typename krylov_traits<Derived>::value_type     value_type;

    protected:
        size_t m_system_size;
        
    public:

        KrylovBase() : m_system_size(derived().system_size()) {}

        /**
         * @brief This function returns an instance of the derived type.
         *
         * @return Derived&
         **/
        Derived &derived()
        {
            return *static_cast<Derived*> ( this );
        }
        
        /**
        * @brief Orthonormalization method.  Applies Arnoldi's method
        *     with modified Gram-Schmidt orthogonalization to the columns of v
        * @param k iteration index
        * @return Upon return, H(k+1,k) = norm(v(:,k+1)) and v(:,k+1) will be orthonormal to all previous v's
        **/
        inline value_type arnoldi ( int k )
        {
            value_type normv = norm ( derived().v ( k+1 ) );
            value_type Hip = value_type ( 0 );
            /// Apply Arnoldi's method with modified Gram-Schmidt orthogonalization to v(k+1)

            for ( unsigned int i = 0; i <= k; ++i )
            {
                derived().H ( i,k ) = dot ( derived().v ( i ),derived().v ( k+1 ) );
                for ( int j = 0; j < m_system_size; ++j )
                    derived().v ( k+1 ) [j] -= derived().H ( i,k ) *derived().v ( i ) [j];
            }
            Hip = norm ( derived().v ( k+1 ) );
            /// Re-orthogonalize if necessary
            if ( normv + value_type ( .001 ) *Hip == normv )
            {
                for ( unsigned int i = 0, nzeros = 0; i <= k; ++i, nzeros+=i )
                {
                    value_type hr = dot ( derived().v ( i ),derived().v ( k+1 ) );
                    derived().H ( i,k ) += hr;
                    for ( int j = 0; j < m_system_size; ++j )
                        derived().v ( k+1 ) [j] -= hr*derived().v ( i ) [j];
                }

                Hip = norm ( derived().v ( k+1 ) );
            }
            /// Normalize
            if ( Hip != value_type ( 0 ) )
                for ( int j = 0; j < m_system_size; ++j )
                    derived().v ( k+1 ) [j] /= Hip;

            return Hip;
        }
        /**
         * \brief Returns the values of the rotation matrix that perform a Givens' rotation.
         *   That is, it returns the following transformation
         *   1/sqrt(dx^2+dy^2)[dx dy] [dx] = [dx^2+dy^2]/sqrt(dx^2+dy^2)
         *   1/sqrt(dx^2+dy^2)[-dy dx][dy]   [0]
         *
         *          /|
         *        /  | dy
         *      /____|
         *        dx
         *
         * \param dx The length of triangle in x.  cos(theta) = dx/sqrt(dx^2+dy^2)
         * \param dy This value will transform to zero. sin(theta) = dy/sqrt(dx^2+dy^2)
         * \param cs dx/sqrt(dx^2+dy^2)
         * \param sn dy/sqrt(dx^2+dy^2)
         *
         **/
        inline void get_rotation ( const value_type &dx, const value_type &dy, value_type &cs, value_type &sn )
        {
            if ( dy == 0.0 )
            {
                cs = 1.0;
                sn = 0.0;
            }
            else if ( std::fabs ( dy ) > std::fabs ( dx ) )
            {
                assert ( dy != 0. );
                value_type temp = dx / dy;
                sn = 1.0 / std::sqrt ( 1.0 + temp * temp );
                cs = temp * sn;
            }
            else
            {
                assert ( dx != 0. );
                value_type temp = dy / dx;
                cs = 1.0 / std::sqrt ( 1.0 + temp * temp );
                sn = temp * cs;
            }
        }

        /**
         * \brief Apply the Givens rotation.  Does the rotation matrix-vector product.
         *
         * \param dx Result of the rotation
         * \param dy Result of the rotation, should equal 0.
         * \param cs Rotation cosine parameter
         * \param sn Rotation sine parameter
         *
         **/
        inline void apply_rotation ( value_type &dx, value_type &dy, const value_type &cs, const value_type &sn )
        {
            value_type temp  =  cs * dx + sn * dy;
            dy = cs * dy - sn * dx;
            dx = temp;
        }

        /**
         * \brief Back substitution solver for an upper triangular linear system. Uy=b. 
         *  Where U = H(0:k,0:k)  and b = g(0:k).
         *
         * \param k Index in the GMRES loop
         * \return Upon return, g(0:k) contains the solution.
         *
         **/
        inline void back_solve ( int k )
        {
            for ( int i = k-1; i >= 0; i-- )
            {
                derived().g ( i ) /= derived().H ( i,i );
                for ( int j = i - 1; j >= 0; j-- )
                    derived().g ( j ) -= derived().H ( j,i ) * derived().g ( i );
            }
        }
        
        inline value_type dot ( const value_type *x, const value_type *y )
        {
            value_type s = value_type ( 0 );
            for ( size_t i = 0; i < m_system_size; ++i )
                s += x[i]*y[i];
            return s;
        }
        
        inline value_type norm ( const value_type *x )
        {
            return std::sqrt ( dot ( x,x ) );
        }
};





#endif
