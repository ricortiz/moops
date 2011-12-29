#ifndef SEMI_IMPLICIT_SDC_RAW_HPP
#define SEMI_IMPLICIT_SDC_RAW_HPP
/**
 * \brief C++ Class.
 * Author: Ricardo Ortiz <ortiz@unc.edu>
 * Carolina Center for Interdisciplinary Applied Mathematics
 * NSF Focused Research Group - Advanced Algorithms and Software for Problems in Computational Bio-Fluid Dynamics
 * $Id:
**/

#include "sdc_base.hpp"
#include "math/ode_solver/euler/backward_euler.hpp"

template<typename function_type>
struct implicit_function_type : public function_type
{
    function_type &m_F;
    implicit_function_type ( function_type &F ) : m_F ( F ) {}

    template<typename value_type>
    inline void operator() ( value_type t, value_type *x, value_type Fx )
    {
        this->Implicit ( t, x, Fx );
    }
};

/**
 * \brief This class implements a fully explicit SDC method.
 * \param implicit_function_type The right hand side function of the differential equation.
 * \param explicit_function_type The right hand side function of the differential equation.
 * \param integrator_type The integrator method used in the correction step.
 * \param sdc_corrections Number of corrections to do.
 **/
template<typename value_type, typename function_type, typename integrator_type, int sdc_nodes, int sdc_corrections>
class SemiImplicitSDC : public SDCBase<SemiImplicitSDC<value_type, function_type, integrator_type, sdc_nodes, sdc_corrections> >
{
    protected:
        typedef BackwardEuler<value_type, function_type> backward_euler_type;

    protected:
        sdc_storage<value_type, sdc_nodes, 0, SDC::SEMI_IMPLICIT> m_storage;
        integrator_type m_integrator;
        size_t m_ode_size;
        function_type   &m_F;
        backward_euler_type m_backward_euler;

    public:
        SemiImplicitSDC ( function_type &F )
            :
            m_storage ( F.ode_size() ),
            m_integrator ( F.ode_size() ),
            m_ode_size ( F.ode_size() ),
            m_F ( F ),
            m_backward_euler ( implicit_function_type<function_type> ( F ) ) {}

        inline const value_type* Fi ( int i ) const { return m_storage.Fi() [i]; }
        inline value_type* Fi ( int i )             { return m_storage.Fi() [i]; }
        inline const value_type** Fi() const         { return m_storage.Fi(); }
        inline value_type** Fi()                     { return m_storage.Fi(); }
        inline const value_type* Fe ( int i ) const { return m_storage.Fe() [i]; }
        inline value_type* Fe ( int i )             { return m_storage.Fe() [i]; }
        inline const value_type** Fe() const         { return m_storage.Fe(); }
        inline value_type** Fe()                     { return m_storage.Fe(); }
        inline const value_type* X ( int i ) const  { return m_storage.X() [i]; }
        inline value_type* X ( int i )              { return m_storage.X() [i]; }
        inline const value_type** X() const          { return m_storage.X(); }
        inline value_type** X()                      { return m_storage.X(); }
        inline void update()                         { m_storage.update(); }
        inline void integrate()                      { m_integrator.integrate ( m_storage.Fi(), m_storage.Fe() ); }
        inline const value_type dt ( int i ) const   { return m_integrator.dt[i]; }
        inline value_type &Immk ( int i, int j ) { return m_integrator.Immk[i][j]; }
        inline const value_type &Immk ( int i, int j ) const { return m_integrator.Immk[i][j]; }
        inline size_t ode_size() { return m_ode_size; }

        inline void init ( value_type *x, value_type *F_i, value_type *F_e ) { m_storage.init ( x, F_i, F_e ); }

        inline void init ( value_type t, value_type *x )
        {
            value_type *x0 = X ( 0 );
            std::copy ( x, x + m_ode_size, x0 );
            m_F.Explicit ( t, x0, Fe ( 0 ) );
            m_F.Implicit ( t, x0, Fi ( 0 ) );
        }

        /**
        * \brief The predictor steps updates xnew and and Fnew by applying forward Euler's method
        *
        * \param xold This is x_k
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
            std::copy ( X ( k ), X ( k ) + m_ode_size, X ( k + 1 ) );
            m_backward_euler ( t, X ( k + 1 ), X ( k ), Fi ( k + 1 ), Fe ( k ), dt );
            m_F.Explicit ( t, X ( k + 1 ), Fe ( k + 1 ) );
        }

        /**
         * \brief The predictor steps updates xnew and and Fnew by applying forward Euler's method
         *
         * \param xold This is x_k
         * \param Fold This is F(x_k)
         * \param xnew The value of x_{k+1}
         * \param Fnew The value of F(x_{k+1})
         * \param t Local time
         * \param dt Local timestep
         * \return void
         *
         * \sa predictor(), corrector()
         **/
        inline void corrector_predictor_step ( const int k, value_type *fdiff, value_type &t, const value_type &dt )
        {
//             fdiff -= Fi ( k+1 );
            std::transform ( fdiff, fdiff + m_ode_size, Fi ( k + 1 ), fdiff, std::minus<value_type>() );
            std::vector<value_type> Fold ( m_ode_size );
            std::copy ( Fe ( k + 1 ), Fe ( k + 1 ) + m_ode_size, Fold.begin() );
            t += dt;
            m_backward_euler ( t, X ( k + 1 ), X ( k ), Fi ( k + 1 ), fdiff, dt );
            m_F.Explicit ( t, X ( k + 1 ), Fe ( k + 1 ) );
            std::transform ( Fe ( k + 1 ), Fe ( k + 1 ) + m_ode_size, Fold.begin(), fdiff, std::minus<value_type>() );
//             fdiff = Fe ( k+1 ) - Fold;
        }



};

template<typename _value_type, typename _function_type, typename _integrator_type, int _sdc_nodes, int _sdc_corrections>
struct sdc_traits<SemiImplicitSDC<_value_type, _function_type, _integrator_type, _sdc_nodes, _sdc_corrections> >
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


