#ifndef EXPLICIT_SDC_RAW_HPP
#define EXPLICIT_SDC_RAW_HPP
/**
 * \brief C++ Class.
 * Author: Ricardo Ortiz <ortiz@unc.edu>
 * Carolina Center for Interdisciplinary Applied Mathematics
 * NSF Focused Research Group - Advanced Algorithms and Software for Problems in Computational Bio-Fluid Dynamics
 * $Id:
**/

#include "sdc_base.hpp"
#include "math/ode_solver/euler/forward_euler.hpp"

/**
 * \brief This class implements a fully explicit SDC method.
 *
 * \param function_type The right hand function of the differential equation.
 * \param integrator_type The integrator method used in the correction step.
 * \param sdc_nodes The umber of subnodes.
 * \param sdc_corrections Number of corrections to do.
 *
 **/
template<typename value_type, typename function_type, typename integrator_type, int sdc_nodes, int sdc_corrections>
class ExplicitSDC : public SDCBase<ExplicitSDC<value_type, function_type, integrator_type, sdc_nodes, sdc_corrections> >
{

    protected:
        typedef ForwardEuler<value_type, function_type> forward_euler_type;

    protected:

        sdc_storage<value_type, sdc_nodes, 0, SDC::EXPLICIT> m_storage;
        integrator_type m_integrator;
        size_t m_ode_size;
        function_type &m_F;
        forward_euler_type m_forward_euler;
        /**< This is the integrator used to compute the spectral integrals of the right hand sides F. **/

    public:

        ExplicitSDC ( function_type &Rhs ) : m_storage ( Rhs.ode_size() ), m_integrator ( Rhs.ode_size() ), m_ode_size ( Rhs.ode_size() ), m_F ( Rhs ), m_forward_euler ( Rhs ) {}

        inline void update()                        { m_storage.update(); }
        inline const value_type *F ( int i ) const  { return m_storage.F() [i]; }
        inline const value_type *X ( int i ) const  { return m_storage.X() [i]; }
        inline value_type *F ( int i )              { return m_storage.F() [i]; }
        inline value_type *X ( int i )              { return m_storage.X() [i]; }
        inline const value_type **F() const         { return m_storage.F(); }
        inline const value_type **X() const         { return m_storage.X(); }
        inline value_type **F()                     { return m_storage.F(); }
        inline value_type **X()                     { return m_storage.X(); }
        inline value_type dt ( int i )              { return m_integrator.dt[i]; }
        inline value_type &Immk ( int i, int j )    { return m_integrator.Immk[i][j]; }
        inline void integrate()                     { m_integrator.integrate ( F() ); }
        inline size_t ode_size()                    { return m_ode_size; }

        inline void init ( const value_type *x, const value_type *Fx ) { m_storage.init ( x, Fx ); }
        inline void init ( value_type t, value_type *x )
        {
            std::copy ( x, x + m_ode_size, X ( 0 ) );
            m_F ( t, X ( 0 ), F ( 0 ) );
        }

        /**
        * \brief The predictor steps updates xnew and and Fnew by applying forward Euler's method
        *
        * \param xold Current iterate x_k
        * \param Fold This is F(x_k)
        * \param xnew The value of x_{k+1}
        * \param Fnew The value of F(x_{k+1})
        * \param t Local substep time
        * \param dt Local substep timestep
        * \return void
        *
        * \sa predictor(), corrector()
        **/
        inline void predictor_step ( const int k, value_type &t, const value_type &dt )
        {
            t += dt;
            m_forward_euler ( X ( k + 1 ), X ( k ), F ( k ), dt );
            m_F ( t, X ( k + 1 ), F ( k + 1 ) );
        }

        /**
         * \brief The predictor steps updates xnew and and Fnew by applying forward Euler's method
         *
         * \param xold This is x_k
         * \param xnew The value of x_{k+1}
         * \param Fnew The value of F(x_{k+1})
         * \param fdiff This is the difference between F(x_{k+1}) and F(x_{k})
         * \param t Local time
         * \param dt Local timestep
         * \return void
         *
         * \sa predictor(), corrector()
         **/
        inline void corrector_predictor_step ( const int k, value_type *fdiff, value_type &t, const value_type &dt )
        {
            value_type Fold[m_ode_size];
            std::copy ( F ( k + 1 ), F ( k + 1 ) + m_ode_size, Fold );
            t += dt;
            m_forward_euler ( X ( k + 1 ), X ( k ), fdiff, dt );
            m_F ( t, X ( k + 1 ), F ( k + 1 ) );
            std::transform ( F ( k + 1 ), F ( k + 1 ) + m_ode_size, Fold, fdiff, std::minus<value_type>() );
        }

};

template<typename _value_type, typename _function_type, typename _integrator_type, int _sdc_nodes, int _sdc_corrections>
struct sdc_traits<ExplicitSDC<_value_type, _function_type, _integrator_type, _sdc_nodes, _sdc_corrections> >
{
    typedef _value_type value_type;
    typedef _integrator_type integrator_type;
    typedef _function_type function_type;
    enum
    {
        sdc_nodes = _sdc_nodes,
        sdc_corrections = _sdc_corrections
    };
};

#endif

