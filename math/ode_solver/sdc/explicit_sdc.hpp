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
 * \param sdc_corrections Number of corrections to do.
 *
 **/
template<typename Derived>
class ExplicitSDC : public SDCBase<ExplicitSDC<Derived> >
{

    protected:
        surface_traits<Derived>::value_type value_type;
        surface_traits<Derived>::function_type function_type;
        surface_traits<Derived>::spectral_intgrator integrator_type;
        enum { sdc_corrections = surface_traits<Derived>::sdc_corrections };
        typedef ForwardEuler<Derived> forward_euler_type;

    private:

        sdc_storage<value_type, integrator_type::sdc_nodes, 0, SDC::EXPLICIT> m_storage;
        integrator_type m_integrator;
        function_type &m_F;
        forward_euler_type m_euler_solver;
        /**< This is the integrator used to compute the spectral integrals of the right hand sides F. **/

    public:

        ExplicitSDC(size_t ode_size) : m_storage(ode_size)
        {
            m_integrator.init(ode_size);
        }

        ExplicitSDC(function_type &rhs) : m_storage(rhs.ode_size()), m_F(rhs), m_euler_solver(rhs)
        {
            m_integrator.init(rhs.ode_size());
        }

        ExplicitSDC(function_type &rhs, size_t ode_size) : m_storage(ode_size), m_F(rhs), m_euler_solver(rhs)
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
        inline size_t ode_size()                 { return m_storage.m_ode_size; }

        inline void init(value_type *x, value_type *Fx) { m_storage.init(x, Fx); }
        inline void setX0(value_type *x) { m_storage.setX0(x); }
        inline void setF0(value_type *Fx) { m_storage.setF0(Fx); }

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
        inline void predictor_step(const int k, value_type &t, const value_type &dt)
        {
            assert(k < integrator_type::sdc_nodes);
            t += dt;
            m_euler_solver(X(k + 1), X(k), F(k), dt);
            m_F(t, X(k + 1), F(k + 1));
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
        inline int corrector_predictor_step(const int k, value_type *fdiff, value_type &t, const value_type &dt)
        {
            assert(k < integrator_type::sdc_nodes);
            std::vector<value_type> Fold(m_storage.m_ode_size, 0.0);
            std::copy(F(k + 1), F(k + 1) + m_storage.m_ode_size, Fold.begin());
            t += dt;
            m_euler_solver(X(k + 1), X(k), fdiff, dt);
            m_F(t, X(k + 1), F(k + 1));
            std::transform(F(k + 1), F(k + 1) + m_storage.m_ode_size, Fold.begin(), fdiff, std::minus<value_type>());
        }

};

template<typename Derived>
struct sdc_traits<ExplicitSDC<Derived> >
{
    surface_traits<Derived>::value_type value_type;
    surface_traits<Derived>::function_type function_type;
    surface_traits<Derived>::spectral_intgrator integrator_type;
    enum
    {
        sdc_nodes = integrator_type::sdc_nodes,
        sdc_corrections = surface_traits<Derived>::sdc_corrections
    };
};

#endif

