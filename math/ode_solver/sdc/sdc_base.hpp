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
        };        

    public:

        SDCBase() {}
        
        /**
         * @brief This function returns an instance of the sdc_method type.
         *
         * @return Derived&
         **/
        Derived &sdc_method()
        {
            return *static_cast<Derived*> ( this );
        }

        /**
         * \brief Predictor loop.
         *
         * \param t Global time
         **/
        inline void predictor (value_type t, value_type Dt )
        {
//             assert ( sdc_method().X() != 0 && sdc_method().F() != 0 && "sdc_base::corrector(): You can not use this method with uninitialized arguments." );
            value_type time = t;

            for ( size_t k = 0; k < m_sdc_nodes - 1; ++k )
            {
                value_type dt = Dt*sdc_method().dt ( k );
                sdc_method().predictor_step ( k, time, dt );
            }
        }

        /**
        * \brief Corrector loop.
        *
        * \param t Global time
        **/
        inline void corrector ( value_type t, value_type Dt )
        {
//             assert ( sdc_method().X() != 0 && sdc_method().F() != 0 && "sdc_base::corrector(): You can not use this method with uninitialized arguments." );
	    size_t ode_size = sdc_method().ode_size();
            std::vector<value_type> fdiff(ode_size,0.0);
            for ( size_t i = 0; i < m_sdc_corrections-1; ++i )
            {
                sdc_method().integrate();
                value_type time = t;
                for ( size_t k = 0; k < m_sdc_nodes - 1; ++k )
                {
                    value_type dt = Dt*sdc_method().dt ( k );
                    for(size_t j = 0; j < ode_size; ++j)
                        fdiff[j] += sdc_method().Immk( k,j ) / sdc_method().dt ( k );
                    sdc_method().corrector_predictor_step ( k, &fdiff[0], time, dt );
                }
                std::fill ( fdiff.begin(),fdiff.end(),value_type ( 0 ) );
            }
            sdc_method().update();
        }

        inline void operator()(value_type t, value_type *x, const value_type *f, value_type dt)
        {
            operator()(t,x,x,f,f,dt);
        }
        
        inline void operator()(value_type t, value_type *x, const value_type *xold, value_type *f, const value_type *fold, value_type dt)
        {
            std::copy(xold,xold+sdc_method().ode_size(),x);
            std::copy(fold,fold+sdc_method().ode_size(),f);
            sdc_method().init(x,f);
            predictor(t,dt);
            corrector(t,dt);
        }

        inline void operator()(value_type t, value_type *x, const value_type *xold, value_type *f1, const value_type *f1old, value_type *f2, const value_type *f2old, value_type dt)
        {
            std::copy(xold,xold+sdc_method().ode_size(),x);
            std::copy(f1old,f1old+sdc_method().ode_size(),f1);
            std::copy(f2old,f2old+sdc_method().ode_size(),f2);
            sdc_method().init(x,f1,f2);
            predictor(t,dt);
            corrector(t,dt);
        }

};


#endif




