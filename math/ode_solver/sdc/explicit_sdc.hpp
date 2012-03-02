#ifndef EXPLICIT_SDC_HPP
#define EXPLICIT_SDC_HPP
/****************************************************************************
** MOOPS -- Modular Object Oriented Particle Simulator
** Copyright (C) 2011-2012  Ricardo Ortiz <ortiz@unc.edu>
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/

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

        internal::sdc_storage<value_type, spectral_integrator_type::sdc_nodes> m_storage;
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

        template<typename function_type>
        inline void operator()(function_type &F, value_type t, value_type *x, value_type *xold, value_type *v, value_type *vold, value_type dt)
        {
            setX0(x);
            setF0(v);
            std::copy(xold, xold + m_storage.ode_size, x);
            std::copy(vold, vold + m_storage.ode_size, v);
            predictor(F,t,dt);
            corrector(F,t,dt);
            update();
        }

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
            m_euler_solver(X(k + 1), X(k), F(k), dt);
            t += dt;
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
        inline void corrector_predictor_step(function_type &G, const int k, value_type *fdiff, value_type &t, const value_type &dt)
        {
            assert(k < spectral_integrator_type::sdc_nodes);
            std::vector<value_type> Fold(m_storage.ode_size, 0.0);
            std::copy(F(k + 1), F(k + 1) + m_storage.ode_size, Fold.begin());
            m_euler_solver(X(k + 1), X(k), fdiff, dt);
            t += dt;
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

