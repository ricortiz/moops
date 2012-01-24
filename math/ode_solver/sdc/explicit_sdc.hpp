#ifndef EXPLICIT_SDC_RAW_HPP
#define EXPLICIT_SDC_RAW_HPP
/**
 * \brief C++ Class.
 * Author: Ricardo Ortiz <ortiz@unc.edu>
 * Carolina Center for Interdisciplinary Applied Mathematics
 * NSF Focused Research Group - Advanced Algorithms and Software for Problems in Computational Bio-Fluid Dynamics
 * $Id:
**/

#include "math/ode_solver/sdc/sdc_base.hpp"
#include "math/ode_solver/sdc/sdc_storage.hpp"
#include "math/ode_solver/euler/forward_euler.hpp"
#include "math/ode_solver/sdc/integrator/spectral_integrator.hpp"

/**
 * \brief This class implements a fully explicit SDC method.
 *
 * \param function_type The right hand function of the differential equation.
 * \param integrator_type The integrator method used in the correction step.
 * \param sdc_corrections Number of corrections to do.
 *
 **/
template<typename value_type, typename spectral_integrator_type = Integrator<value_type,gauss_lobatto>, int sdc_corrections = 8>
class ExplicitSDC : public SDCBase<ExplicitSDC<value_type,spectral_integrator_type,sdc_corrections> >
{

    private:

        sdc_storage<value_type, spectral_integrator_type::sdc_nodes, 0, SDC::EXPLICIT> m_storage;
        spectral_integrator_type m_integrator;
        ForwardEuler m_euler_solver;
        /**< This is the integrator used to compute the spectral integrals of the right hand sides F. **/

    public:

        ExplicitSDC(size_t ode_size) : m_storage(ode_size), m_euler_solver(ode_size)
        {
            m_integrator.init(ode_size);
        }

        
        inline void update()                     { m_storage.update(); }
        inline const value_type *F(int i) const  { return m_storage.F()[i]; }
        inline const value_type *X(int i) const  { return m_storage.X()[i]; }
        inline value_type *F(int i)              { return m_storage.F()[i]; }
        inline value_type *X(int i)              { return m_storage.X()[i]; }
        inline const value_type **F() const      { return m_storage.F(); }
        inline const value_type **X() const      { return m_storage.X(); }
        inline value_type **F()                  { return m_storage.F(); }
        inline value_type **X()                  { return m_storage.X(); }
        inline value_type dt(int i)              { return m_integrator.dt(i); }
        inline value_type &Immk(int i, int j)    { return m_integrator.Immk[i][j]; }
        inline void integrate(value_type Dt)     { m_integrator.integrate(F(), Dt); }
        inline size_t ode_size()                 { return m_storage.ode_size; }

        inline void init(value_type *x, value_type *Fx) { m_storage.init(x, Fx); }
        inline void setX0(value_type *x) { m_storage.setX0(x); }
        inline void setF0(value_type *Fx) { m_storage.setF0(Fx); }

//         (function_type &F, value_type t, value_type *x, value_type *v, value_type dt)

        template<typename function_type>
        inline void operator()(function_type &F, value_type t, value_type *x, value_type *v, value_type dt)
        {
            setX0(x);
            setF0(v);
            predictor(F,t,dt);
            corrector(F,t,dt);
            update();
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
        template<typename function_type>
        inline void predictor_step(function_type &G, const int k, value_type &t, const value_type &dt)
        {
            assert(k < spectral_integrator_type::sdc_nodes);
            t += dt;
            m_euler_solver(X(k + 1), X(k), F(k), dt);
            G(t, X(k + 1), F(k + 1));
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
        template<typename function_type>
        inline int corrector_predictor_step(function_type &G, const int k, value_type *fdiff, value_type &t, const value_type &dt)
        {
            assert(k < spectral_integrator_type::sdc_nodes);
            std::vector<value_type> Fold(m_storage.ode_size, 0.0);
            std::copy(F(k + 1), F(k + 1) + m_storage.ode_size, Fold.begin());
            t += dt;
            m_euler_solver(X(k + 1), X(k), fdiff, dt);
            G(t, X(k + 1), F(k + 1));
            std::transform(F(k + 1), F(k + 1) + m_storage.ode_size, Fold.begin(), fdiff, std::minus<value_type>());
        }

};

template<typename _value_type, typename _spectral_integrator_type, int _sdc_corrections>
struct sdc_traits<ExplicitSDC<_value_type,_spectral_integrator_type,_sdc_corrections> >
{
    typedef _value_type value_type;
    enum
    {
        sdc_nodes = _spectral_integrator_type::sdc_nodes,
        sdc_corrections = _sdc_corrections
    };
};

#endif

