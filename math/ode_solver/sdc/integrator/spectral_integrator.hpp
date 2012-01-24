#ifndef SPECTRAL_INTEGRATOR_HPP
#define SPECTRAL_INTEGRATOR_HPP

#include<vector>
#include "quadrature_rule.hpp"

template<typename value_type, quadrature_type quadrature, int _sdc_nodes, int _multirate_nodes, int _total_nodes>
struct SpectralIntegrator
        : public QuadratureRule<value_type, quadrature, _sdc_nodes, _multirate_nodes>::type
{
    enum
    {
        sdc_nodes = _sdc_nodes,
        multirate_nodes = _multirate_nodes,
        total_nodes = _total_nodes
    };
    value_type &Immk;
    SpectralIntegrator<value_type, quadrature, total_nodes, 2, total_nodes> sdc_integrator;

    SpectralIntegrator() : Immk(sdc_integrator.Immk)
    {
    }

    void init(size_t _ode_size)
    {
        sdc_integrator.init(_ode_size);
    }

    inline void integrate(const value_type *F1, const value_type *F2, value_type Dt)
    {
        sdc_integrator.integrate(F1, Dt);
        for(int i = 0; i < total_nodes - 1; ++i)
            for(int j = 0; j < sdc_nodes; ++j)
            {
                value_type c = this->matrix(i, j) * Dt;
                for(size_t k = 0; k < sdc_integrator.ode_size; ++k)
                    Immk[i][k] += F2[j][k] * c;
            }
    }
};


template<typename value_type, quadrature_type quadrature, int _sdc_nodes>
struct SpectralIntegrator<value_type, quadrature, _sdc_nodes, 2, _sdc_nodes>
        : public QuadratureRule<value_type, quadrature, _sdc_nodes, 0>::type
{
    enum
    {
        sdc_nodes = _sdc_nodes
    };

    std::vector<value_type> Immk[sdc_nodes - 1];
    size_t ode_size;

    void init(size_t _ode_size)
    {
        for(int i = 0; i < sdc_nodes-1; ++i)
            Immk[i].resize(_ode_size,0.0);
        ode_size = _ode_size;
    }

    inline void integrate(value_type **F, value_type Dt)
    {
        for(int i = 0; i < sdc_nodes - 1; ++i)
        {
            std::fill(Immk[i].begin(), Immk[i].end(), value_type(0));
            for(int j = 0; j < sdc_nodes; ++j)
            {
                value_type c = this->matrix(i,j) * Dt;
                for(size_t k = 0; k < ode_size; ++k)
                    Immk[i][k] += F[j][k] * c;
            }
        }
    }

    inline void integrate(value_type **F1, value_type **F2, value_type Dt)
    {
        // TODO: Fix me, inefficient!
        value_type *F[sdc_nodes];
        for(int i = 0; i < sdc_nodes; ++i)
            F[i] = new value_type[ode_size];
        for(int i = 0; i < sdc_nodes; ++i)
        {
            for(size_t k = 0; k < ode_size; ++k)
                F[i][k] = F1[i][k] + F2[i][k];
            //             std::copy(F[i],F[i]+ode_size,std::ostream_iterator<value_type>(std::cout, " "));std::cout << std::endl;
        }
        integrate(F, Dt);
        for(size_t i = 0; i < sdc_nodes; ++i)
            delete [] F[i];
    }
    inline void integrate(const value_type **F1, const value_type **F2, const value_type **F3, value_type Dt)
    {
        // TODO: Fix me, inefficient!
        value_type F[sdc_nodes][ode_size];
        for(int i = 0; i < sdc_nodes; ++i)
            for(size_t k = 0; k < ode_size; ++k)
                F[i][k] = F1[i][k] + F2[i][k];
        integrate(F, F3, Dt);
    }
    inline void integrate(const value_type **F1, const value_type **F2, const value_type **F3, const value_type **F4, value_type Dt)
    {
        // TODO: Fix me, inefficient!
        value_type F[sdc_nodes][ode_size];
        for(int i = 0; i < sdc_nodes; ++i)
            for(size_t k = 0; k < ode_size; ++k)
                F[i][k] = F1[i][k] + F2[i][k];
        integrate(F, F3, F4, Dt);
    }
    inline void integrate(const value_type **F1, const value_type **F2, const value_type **F3, const value_type **F4, const value_type **F5, value_type Dt)
    {
        // TODO: Fix me, inefficient!
        value_type F[sdc_nodes][ode_size];
        for(int i = 0; i < sdc_nodes; ++i)
            for(size_t k = 0; k < ode_size; ++k)
                F[i][k] = F1[i][k] + F2[i][k];
        integrate(F, F3, F4, F5, Dt);
    }
};

template < typename value_type, quadrature_type quadrature = other, int sdc_nodes = 5, int multirate_nodes = 2, int total_nodes = (sdc_nodes - 1) * (multirate_nodes - 1) + 1 >
struct Integrator : public SpectralIntegrator<value_type,quadrature,sdc_nodes,multirate_nodes,total_nodes>
    {};

#endif
