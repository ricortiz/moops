#ifndef SEMI_IMPLICIT_SDC_RAW_HPP
#define SEMI_IMPLICIT_SDC_RAW_HPP
/**
 * \brief C++ Class.
 * Author: Ricardo Ortiz <ortiz@unc.edu>
 * Carolina Center for Interdisciplinary Applied Mathematics
 * NSF Focused Research Group - Advanced Algorithms and Software for Problems in Computational Bio-Fluid Dynamics
 * $Id:
**/

#include "math/ode_solver/sdc/sdc_base.hpp"
#include "math/ode_solver/sdc/sdc_storage.hpp"
#include "math/ode_solver/euler/backward_euler.hpp"
#include "math/ode_solver/euler/forward_euler.hpp"
#include "math/ode_solver/sdc/integrator/spectral_integrator.hpp"

template<typename function_type>
struct implicit_function
{
    function_type &m_F;
    implicit_function(function_type &F) : m_F(F) {}

    template<typename value_type>
    inline void operator()(value_type t, value_type *x, value_type *Fx)
    {
        m_F.Implicit(t, x, Fx);
    }
};

/**
 * \brief This class implements a fully explicit SDC method.
 * \param implicit_function_type The right hand side function of the differential equation.
 * \param explicit_function_type The right hand side function of the differential equation.
 * \param integrator_type The integrator method used in the correction step.
 * \param sdc_corrections Number of corrections to do.
 **/
template < typename value_type, typename spectral_integrator_type = Integrator<value_type, gauss_lobatto>, int sdc_corrections = 8 >
class SemiImplicitSDC : public SDCBase<SemiImplicitSDC<value_type, spectral_integrator_type, sdc_corrections> >
{
    protected:
        typedef BackwardEuler<value_type> backward_euler_type;
        typedef ForwardEuler     forward_euler_type;

    protected:
        internal::sdc_storage < value_type,
        spectral_integrator_type::sdc_nodes,
        internal::SEMI_IMPLICIT >                    m_storage;
        spectral_integrator_type                                m_integrator;
        backward_euler_type                                     m_backward_euler;
        forward_euler_type                                      m_forward_euler;

    public:
        SemiImplicitSDC(size_t ode_size) : m_storage(ode_size), m_backward_euler(ode_size), m_forward_euler(ode_size)
        {
            m_integrator.init(ode_size);

        }

        inline const value_type* Fi(int i) const { return m_storage.Fi()[i]; }
        inline value_type* Fi(int i)             { return m_storage.Fi()[i]; }
        inline const value_type** Fi() const         { return m_storage.Fi(); }
        inline value_type** Fi()                     { return m_storage.Fi(); }
        inline const value_type* Fe(int i) const { return m_storage.Fe()[i]; }
        inline value_type* Fe(int i)             { return m_storage.Fe()[i]; }
        inline const value_type** Fe() const         { return m_storage.Fe(); }
        inline value_type** Fe()                     { return m_storage.Fe(); }
        inline const value_type* X(int i) const  { return m_storage.X()[i]; }
        inline value_type* X(int i)              { return m_storage.X()[i]; }
        inline const value_type** X() const          { return m_storage.X(); }
        inline value_type** X()                      { return m_storage.X(); }
        inline void update()                         { m_storage.update(); }
        inline void integrate(value_type Dt)                      { m_integrator.integrate(m_storage.Fi(), m_storage.Fe(), Dt); }
        inline value_type dt(int i)              { return m_integrator.dt(i); }
        inline value_type &Immk(int i, int j) { return m_integrator.Immk[i][j]; }
        inline const value_type &Immk(int i, int j) const { return m_integrator.Immk[i][j]; }
        inline size_t ode_size() { return m_storage.ode_size; }

        template<typename function_type>
        inline void operator()(function_type &F, value_type t, value_type *x, value_type *v, value_type dt)
        {
            m_storage.setX0(x);
            predictor(F, t, dt);
            corrector(F, t, dt);

            update();
            std::transform(Fe(0), Fe(0) + m_storage.ode_size, Fi(0), v, std::plus<value_type>());
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
        template<typename function_type>
        inline void predictor_step(function_type &F, const int k, value_type &t, const value_type &dt)
        {
            assert(k < spectral_integrator_type::sdc_nodes);
            implicit_function<function_type> G(F);
            t += dt;
            std::vector<value_type> buffer(m_storage.ode_size);
            m_forward_euler(&buffer[0], X(k), Fe(k), dt);
            m_backward_euler(G, t, X(k + 1), &buffer[0], Fi(k + 1), dt);
            F.Explicit(t, X(k + 1), Fe(k + 1));
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
        template<typename function_type>
        inline void corrector_predictor_step(function_type &F, const int k, value_type *fdiff, value_type &t, const value_type &dt)
        {
            assert(k < spectral_integrator_type::sdc_nodes);
            implicit_function<function_type> G(F);
            std::vector<value_type> buffer(m_storage.ode_size);
	    
            t += dt;
            std::transform(fdiff, fdiff + m_storage.ode_size, Fi(k + 1), fdiff, std::minus<value_type>());
            m_forward_euler(&buffer[0], X(k), fdiff, dt);
            m_backward_euler(G, t, X(k + 1), &buffer[0], Fi(k + 1), dt);
	    
            std::copy(Fe(k + 1), Fe(k + 1) + m_storage.ode_size, buffer.begin());
            F.Explicit(t, X(k + 1), Fe(k + 1));
            std::transform(Fe(k + 1), Fe(k + 1) + m_storage.ode_size, buffer.begin(), fdiff, std::minus<value_type>());
        }



};

template<typename _value_type, typename _integrator_type, int _sdc_corrections>
struct sdc_traits<SemiImplicitSDC<_value_type, _integrator_type, _sdc_corrections> >
{
    typedef _value_type value_type;
    typedef _integrator_type integrator_type;
    enum
    {
        sdc_nodes = _integrator_type::sdc_nodes,
        sdc_corrections = _sdc_corrections
    };
};

#endif


