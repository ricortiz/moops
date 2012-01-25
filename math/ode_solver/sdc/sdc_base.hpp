#ifndef SDC_BASE_RAW_HPP
#define SDC_BASE_RAW_HPP
/**
 * @brief C++ Class.
 * Author: Ricardo Ortiz <ortiz@unc.edu>
 * Carolina Center for Interdisciplinary Applied Mathematics
 * NSF Focused Research Group - Advanced Algorithms and Software for Problems in Computational Bio-Fluid Dynamics
 * $Id:
**/
#include<numeric>
#include<cassert>
#include<stdexcept>


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
            sdc_nodes = sdc_traits<Derived>::sdc_nodes,
            /**< The number of SDC nodes at compile-time. This is just a copy of the value provided
            * by the \a Derived type. **/

            sdc_corrections = sdc_traits<Derived>::sdc_corrections
            /**< The number of SDC corrections sweps at compile-time. This is just a copy of the value provided
            * by the \a Derived type. **/
        };

        value_type m_residuals[sdc_corrections][sdc_nodes - 1];

    public:

        /**
         * @brief This function returns an instance of the sdc_method type.
         *
         * @return Derived&
         **/
        Derived &sdc_method()
        {
            return *static_cast<Derived *>(this);
        }

        /**
         * \brief Predictor loop.
         *
         * \param t Global time
         **/
        template<typename function_type>
        inline void predictor(function_type &F, value_type t, value_type Dt)
        {
//             assert ( sdc_method().X() != 0 && sdc_method().F() != 0 && "sdc_base::corrector(): You can not use this method with uninitialized arguments." );
            value_type time = t;

            for(size_t k = 0; k < sdc_nodes - 1; ++k)
            {
                value_type dt = Dt * sdc_method().dt(k);
                sdc_method().predictor_step(F, k, time, dt);
                check_convergence(0, k);
            }
        }

        /**
        * \brief Corrector loop.
        *
        * \param t Global time
        **/
        template<typename function_type>
        inline void corrector(function_type &F, value_type t, value_type Dt)
        {
//             assert ( sdc_method().X() != 0 && sdc_method().F() != 0 && "sdc_base::corrector(): You can not use this method with uninitialized arguments." );
            size_t ode_size = sdc_method().ode_size();
            std::vector<value_type> fdiff(ode_size, 0.0);

            for(size_t i = 0; i < sdc_corrections - 1; ++i)
            {
                sdc_method().integrate(Dt);
                value_type time = t;
                for(size_t k = 0; k < sdc_nodes - 1; ++k)
                {
                    value_type dt = Dt * sdc_method().dt(k);
                    for(size_t j = 0; j < ode_size; ++j)
                        fdiff[j] += sdc_method().Immk(k, j) / dt;
                    sdc_method().corrector_predictor_step(F, k, &fdiff[0], time, dt);
                    check_convergence(i + 1, k);
                }
//                 if(m_residuals[i][sdc_nodes - 2] < 1e-13)
//                 {
//                     sdc_method().update();
//                     return;
//                 }
                std::fill(fdiff.begin(), fdiff.end(), value_type(0));
            }
//             sdc_method().update();
            for ( size_t i = 0; i < sdc_corrections; ++i )
            {
                for ( size_t k = 0; k < sdc_nodes - 1; ++k )
                    std::cout << m_residuals[i][k] << " ";
                std::cout << std::endl;
            }
            std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        }

        void check_convergence(int i, int k)
        {
            size_t ode_size = sdc_method().ode_size();
            std::vector<value_type> tmp(ode_size, 0.0);
            std::transform(sdc_method().X(k), sdc_method().X(k) + ode_size, &sdc_method().Immk(k, 0), tmp.begin(), std::plus<value_type>());
            std::transform(sdc_method().X(k + 1), sdc_method().X(k + 1) + ode_size, tmp.begin(), tmp.begin(), std::minus<value_type>());
            m_residuals[i][k] = std::sqrt(std::inner_product(tmp.begin(), tmp.end(), tmp.begin(), 0.0));
            if(m_residuals[i][k] > 5)
            {
                std::cout << "correction = " << i << ", iteration = " << k << std::endl;
                for(int j = 0; j <= i; ++j)
                {
                    for(int l = 0; l <= k; ++l)
                        std::cout << m_residuals[j][l] << " ";
                    std::cout << std::endl;
                }
                throw std::out_of_range("sdc residual is out of bound: SDC_BASE::CHECK_CONVERGENCE");
            }
        }

};


#endif




