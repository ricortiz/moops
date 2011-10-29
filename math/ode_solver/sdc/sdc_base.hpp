#ifndef SDC_BASE_RAW_HPP
#define SDC_BASE_RAW_HPP
/**
 * @brief C++ Class.
 * Author: Ricardo Ortiz <ortiz@unc.edu>
 * Carolina Center for Interdisciplinary Applied Mathematics
 * NSF Focused Research Group - Advanced Algorithms and Software for Problems in Computational Bio-Fluid Dynamics
 * $Id:
**/

#include<cassert>
#include "math/nonlinear_solver/inexact_newton.hpp"
#include "sdc_storage.hpp"

template<typename T> struct sdc_traits;

/**
 * \brief SDC base class.  Base class for all SDC variants.
 *
 * This class is the base that is inherited by all SDC types.
 * Most of the SDC API is contained in this class.
 *
 * \param Derived is the sdc_method type, e.g. an explicit sdc type, etc.
 *

 * When writing a function taking SDC objects as argument, if you want your function
 * to take as argument any sdc_method, just let it take a SDCBase argument.
 * As an example, here is a function time_stepping which, given
 * an SDC subclass, computes an iteration of the SDC method.
 *
 * \code
   template<typename Derived>
   void time_stepping(const SDCBase<Derived>& sdc)
   {
     ...
     sdc.predictor(x,F,t);
     sdc.corrector(x,F,t);
   }
 * \endcode
 * The Derived type should define the following two methods used in this class:
 * predictor_step(), corrector_step()
 **/
template<typename Derived>
class SDCBase
{

    public:
        typedef Derived sdc_method_type;
        typedef typename sdc_traits<Derived>::value_type value_type;
    private:

        enum
        {
            m_sdc_nodes = sdc_traits<Derived>::sdc_nodes,
            /**< The number of SDC nodes at compile-time. This is just a copy of the value provided
            * by the \a Derived type. **/

            m_sdc_corrections = sdc_traits<Derived>::sdc_corrections,
            /**< The number of SDC corrections sweps at compile-time. This is just a copy of the value provided
            * by the \a Derived type. **/

            m_ode_size = sdc_traits<Derived>::ode_size
            /**< Size of the problem at compile-time. This is just a copy of the value provided
             * by the \a Derived type. **/
        };


    public:

        /**
         * @brief This function returns an instance of the sdc_method type.
         *
         * @return Derived&
         **/
        Derived &sdc_method()
        {
            return *static_cast<Derived*> ( this );
        }
        
//         inline const value_type *F ( int i ) const  { return sdc_method().F(i); }
//         inline const value_type *X ( int i ) const  { return sdc_method().X(i); }
//         inline value_type *F ( int i )              { return sdc_method().F(i); }
//         inline value_type *X ( int i )              { return sdc_method().X(i); }

        /**
         * \brief Predictor loop.
         *
         * \param x Vector containing step's data for unknown variables
         * \param F Vector containing step's data for rigth hand side
         * \param t Global time
         **/
        template<typename function_type>
        inline void predictor (function_type &V, value_type t, value_type Dt )
        {
//             assert ( sdc_method().X() != 0 && sdc_method().F() != 0 && "sdc_base::corrector(): You can not use this method with uninitialized arguments." );
            value_type time = t;

            for ( size_t k = 0; k < m_sdc_nodes - 1; ++k )
            {
                value_type dt = Dt*sdc_method().dt ( k );
                sdc_method().predictor_step ( V, k, time, dt );
            }
        }

        /**
        * \brief Corrector loop.
        *
        * \param x Vector containing step's data for unknown variables
        * \param F Vector containing step's data for rigth hand side
        * \param t Global time
        **/
        template<typename function_type>
        inline void corrector ( function_type&V, value_type t, value_type Dt )
        {
//             assert ( sdc_method().X() != 0 && sdc_method().F() != 0 && "sdc_base::corrector(): You can not use this method with uninitialized arguments." );
            value_type fdiff[m_ode_size];
            for ( size_t i = 0; i < m_sdc_corrections; ++i )
            {
                std::fill ( fdiff,fdiff+m_ode_size,value_type ( 0 ) );
                sdc_method().integrate();
                value_type time = t;
                for ( size_t k = 0; k < m_sdc_nodes - 1; ++k )
                {
                    value_type dt = Dt*sdc_method().dt ( k );
                    for(size_t j = 0; j < m_ode_size; ++j)
                        fdiff[j] += sdc_method().Immk( k,j ) / sdc_method().dt ( k );
                    sdc_method().corrector_predictor_step ( V, k, fdiff, time, dt );
                }
            }

            sdc_method().update();
        }

        /**
        * \brief Forward Euler step, used in predictor and corrector steps.
        *
        * \param xnew New update
        * \param xold Old value
        * \param F Function value: F(x_k)
        * \param dt local timestep
        *
        * \sa predictor_step()
        **/
        inline void forward_euler ( value_type *xnew, const value_type *xold, const value_type *F, const value_type &dt )
        {
            for ( size_t i = 0; i < m_ode_size; ++i )
                xnew[i] = xold[i] + dt * F[i];
        }

        /**
         * \brief Backward Euler step, used in predictor and corrector steps.
         *
         * \param xnew New update
         * \param xold Old value
         * \param F Function value: F(x_k)
         * \param dt local timestep
         *
         * \sa predictor_step()
         **/
        template<typename operator_type>
        inline void backward_euler ( value_type *x, value_type t, const value_type *xold, const value_type *Fold, value_type *Fi, value_type dt, operator_type &Vi )
        {
            InexactNewtonMethod<value_type,m_ode_size> newton_solve;
            value_type rhs[m_ode_size];
            std::fill(rhs,rhs+m_ode_size,value_type(0));
            
            forward_euler ( rhs, xold, Fold, dt );
            implicit_operator<operator_type,value_type> F ( Vi,t,dt,rhs );
            newton_solve ( F, x, 1e-12, 1e-12 );
            for ( size_t i = 0; i < m_ode_size; ++i )
                Fi[i] = ( x[i] - rhs[i] ) / dt;
        }

        template<typename function_type, typename value_type>
        struct implicit_operator
        {
            function_type &m_V;
            const value_type &m_t;
            const value_type &m_dt;
            const value_type *m_rhs;

            implicit_operator ( function_type &V, const value_type &t, const value_type &dt, const value_type *rhs )
                    : m_V ( V ), m_t ( t ), m_dt ( dt ), m_rhs ( rhs ) {}

            void operator() ( const value_type *x, value_type *Fx )
            {
                value_type Vx[m_ode_size];
                m_V.I ( m_t, x, Vx );
                for ( size_t i = 0; i < m_ode_size; ++i )
                    Fx[i] = x[i] - m_dt * Vx[i] - m_rhs[i];
            }

        };
};


#endif




