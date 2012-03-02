#ifndef GENERALIZED_MINIMAL_RESIDUAL_METHOD_HPP
#define GENERALIZED_MINIMAL_RESIDUAL_METHOD_HPP
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
#include<algorithm>

#include "krylov_storage.hpp"
#include "krylov_base.hpp"

template<typename value_type, int krylov_space_maximum_dimension>
class GeneralizedMinimalResidualMethod : public KrylovBase<GeneralizedMinimalResidualMethod<value_type, krylov_space_maximum_dimension> >
{
    protected:
        krylov_storage<value_type, krylov_space_maximum_dimension> m_storage;

    private:
        size_t m_system_size;
        enum
        {
            k_max = krylov_space_maximum_dimension
        };

    public:

        GeneralizedMinimalResidualMethod(size_t system_size) : m_storage(system_size), m_system_size(system_size){}
        inline value_type &H(size_t i, size_t j) { return m_storage.H(i, j); }
        inline value_type *H() { return m_storage.H(); }
        inline value_type *residual() { return m_storage.residual(); }
        inline value_type &residual(size_t i) { return m_storage.residual(i); }
        inline value_type &s(size_t i) { return m_storage.s(i); }
        inline value_type &c(size_t i) { return m_storage.c(i); }
        inline value_type &g(size_t i) { return m_storage.g(i); }
        inline value_type *g() { return m_storage.g(); }
        inline value_type *v(size_t i) { return m_storage.v(i); }
        inline size_t system_size() { return m_system_size; }

        /**
        * \brief Gneralized Minimal Residual Method.  Solves the linear system:  F(x)= b, where F(x) = A*x.
        *
        * \param F This is a linear operator on x.
        * \param b Right hand side vector of the linear system.
        * \param x Initial guess and solution vector
        * \param stats This is an optional vector to collect statistics of the method.  Defaults to 0.
        * \return norm of the residual
        **/
        template<typename operator_type, typename vector_type>
        inline value_type operator()(operator_type &F, const value_type *b, value_type *x, value_type errtol = 1e-6, std::map<std::string, vector_type> *stats = 0)
        {
            /// Compute the residual
            F(x, residual());
            std::transform(b, b + m_system_size, residual(), residual(), std::minus<value_type>());
            value_type rho = this->norm(residual());
            errtol *= this->norm(b);
            if(stats)
            {
                std::map<std::string, vector_type> &s = *stats;
                std::string residuals_gmres_key = "gmres_residuals";
                s[residuals_gmres_key].push_back(rho);
            }

            /// Test for early termination
            if(rho < errtol)
            {
                if(stats)
                {
                    std::map<std::string, vector_type> &s = *stats;
                    std::string fn_eval_gmres_key = "gmres_fn_eval";
                    s[fn_eval_gmres_key].push_back(1);
                }
                return rho;
            }
            std::fill(g(), g() + k_max + 1, value_type(0));
            g(0) = rho;

            for(size_t i = 0 ; i < m_system_size; ++i)
                v(0)[i] = residual(i) / rho;
            /// Start GMRES iteration
            int k;
            for(k = 0; k < k_max && rho > errtol; ++k)
            {
                /// Evaluate linear operator
                F(v(k), v(k + 1));

                /// Apply Arnoldi's method with modified Gram-Schmidt orthogonalization to the colums of v
                value_type ld = this->arnoldi(k);

                /// Perform Givens rotations in order to eliminate the H(k+1,k) entry from the
                /// Hessenberg matrix.
                for(int i = 0; i < k; ++i)
                    this->apply_rotation(H(i, k), H(i + 1, k), c(i), s(i));
                this->get_rotation(H(k, k), ld, c(k), s(k));
                this->apply_rotation(H(k, k), ld, c(k), s(k));
                this->apply_rotation(g(k), g(k + 1), c(k), s(k));

                assert(ld < 1e-10 || !"rotation failed: Lower diagonal is non-zero");

                ///  The norm of the residual is the last entry in g
                rho = std::fabs(g(k + 1));
                if(stats)
                {
                    std::map<std::string, vector_type> &s = *stats;
                    std::string residuals_gmres_key = "gmres_residuals";
                    s[residuals_gmres_key].push_back(rho);
                }
            }
            if(stats)
            {
                std::map<std::string, vector_type> &s = *stats;
                std::string fn_eval_gmres_key = "gmres_fn_eval";
                s[fn_eval_gmres_key].push_back(k + 1);
            }
            
            /// Solve the least square problem by solving the upper triangular linear system
            this->back_solve(k);
            
            /// Update the solution
            for(int i = 0; i < k; i++)
            {
                for(size_t j = 0; j < m_system_size; ++j)
                    x[j] += g(i)*v(i)[j];
            }

            return rho;
        }
};

template<typename _value_type, size_t _krylov_space_max_dim>
struct krylov_traits<GeneralizedMinimalResidualMethod<_value_type, _krylov_space_max_dim> >
{
    typedef _value_type value_type;
    enum
    {
        krylov_space_max_dim = _krylov_space_max_dim
    };
};

#endif
