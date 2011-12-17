#ifndef CLENSHAW_CURTIS_RAW_HPP
#define CLENSHAW_CURTIS_RAW_HPP

#include "clenshaw_curtis_coefficients.hpp"

template<typename value_type, int _ode_size, int _sdc_nodes, int _multirate_nodes, int _total_nodes>
struct SDCSpectralIntegrator;

template<typename value_type, size_t _ode_size, size_t _sdc_nodes, size_t _multirate_nodes, size_t _total_nodes>
struct SDCSpectralIntegrator<value_type, _ode_size,_sdc_nodes,_multirate_nodes,_total_nodes>
            : public SdcCoefficients<_sdc_nodes,_multirate_nodes>::type
{

    enum
    {
        ode_size = _ode_size,
        sdc_nodes = _sdc_nodes,
        multirate_nodes = _multirate_nodes,
        total_nodes = _total_nodes
    };

    value_type *Immk;
    SDCSpectralIntegrator<value_type,ode_size,total_nodes,2,total_nodes> sdc_integrator;
    value_type coeff[sdc_nodes* ( total_nodes-1 ) ];
    value_type *points;
    value_type *dt;

    SDCSpectralIntegrator() : Immk ( sdc_integrator.Immk ), points ( sdc_integrator.points ), dt ( sdc_integrator.dt )
    {
        get_coefficients ( coeff );
        value_type h = M_PI / ( sdc_nodes - 1 );
    }

    inline void integrate ( const value_type *F1, const value_type *F2 )
    {
        sdc_integrator.integrate ( F1 );
        for ( size_t j = 0; j < total_nodes - 1; ++j )
            for ( size_t i = 0; i < sdc_nodes; ++i )
                for ( size_t k = 0; k < ode_size; ++k )
                    Immk[j][k] += F2[i][k] * coeff[sdc_nodes*j+i];
    }
};

template<typename value_type, int _ode_size, int _sdc_nodes>
struct SDCSpectralIntegrator<value_type, _ode_size,_sdc_nodes,2,_sdc_nodes> : public SdcCoefficients<_sdc_nodes,0>::type
{

    enum
    {
        ode_size = _ode_size,
        sdc_nodes = _sdc_nodes
    };

    value_type coeff[sdc_nodes* ( sdc_nodes-1 ) ];
    value_type points[sdc_nodes];
    value_type dt[sdc_nodes-1];
    value_type Dt;

    value_type* Immk[sdc_nodes-1];

    SDCSpectralIntegrator()
    {
        for ( int i = 0; i < sdc_nodes-1; ++i )
            Immk[i] = new value_type[ode_size];
        get_coefficients ( coeff );
        value_type h = M_PI / ( sdc_nodes - 1 );
        for ( size_t i = 0; i < sdc_nodes; ++i )
        {
            value_type theta = -M_PI + i * h;
            points[i] = .5 * ( std::cos ( theta ) + 1 );
        }
        for ( size_t i = 0; i < sdc_nodes-1; ++i )
            dt[i] = points[i+1]-points[i];
    }
    ~SDCSpectralIntegrator()
    {
        for ( int i = 0; i < sdc_nodes-1; ++i )
            delete [] Immk[i];
    }
    inline void integrate ( value_type **F )
    {
        for ( size_t j = 0; j < sdc_nodes - 1; ++j )
        {
            std::fill ( Immk[j],Immk[j]+ode_size,value_type ( 0 ) );
            for ( size_t i = 0; i < sdc_nodes; ++i )
            {
//                 std::copy(F[i],F[i]+ode_size,std::ostream_iterator<value_type>(std::cout, " "));
                for ( size_t k = 0; k < ode_size; ++k )
                    Immk[j][k] += F[i][k] * coeff[sdc_nodes*j+i];
            }
        }
    }

    inline void integrate ( value_type **F1, value_type **F2 )
    {
        value_type *F[sdc_nodes];
        for(size_t i = 0; i < sdc_nodes; ++i)
            F[i] = new value_type[ode_size];
        for ( size_t i = 0; i < sdc_nodes; ++i )
        {            
            for ( size_t k = 0; k < ode_size; ++k )
                F[i][k] = F1[i][k]+F2[i][k];
//             std::copy(F[i],F[i]+ode_size,std::ostream_iterator<value_type>(std::cout, " "));std::cout << std::endl;
        }
        integrate ( F );
        for(size_t i = 0; i < sdc_nodes; ++i)
            delete [] F[i];
    }
    inline void integrate ( const value_type **F1, const value_type **F2, const value_type **F3 )
    {
        value_type F[sdc_nodes][ode_size];
        for ( size_t i = 0; i < sdc_nodes; ++i )
            for ( size_t k = 0; k < ode_size; ++k )
                F[i][k] = F1[i][k]+F2[i][k];
        integrate ( F, F3 );
    }
    inline void integrate ( const value_type **F1, const value_type **F2, const value_type **F3, const value_type **F4 )
    {
        value_type F[sdc_nodes][ode_size];
        for ( size_t i = 0; i < sdc_nodes; ++i )
            for ( size_t k = 0; k < ode_size; ++k )
                F[i][k] = F1[i][k]+F2[i][k];
        integrate ( F, F3, F4 );
    }
    inline void integrate ( const value_type **F1, const value_type **F2, const value_type **F3, const value_type **F4, const value_type **F5 )
    {
        value_type F[sdc_nodes][ode_size];
        for ( size_t i = 0; i < sdc_nodes; ++i )
            for ( size_t k = 0; k < ode_size; ++k )
                F[i][k] = F1[i][k]+F2[i][k];
        integrate ( F, F3, F4, F5 );
    }
};


template<typename value_type,  size_t _sdc_nodes, size_t _multirate_nodes, size_t _total_nodes>
struct SDCSpectralIntegrator<value_type, 0,_sdc_nodes,_multirate_nodes,_total_nodes>
: public SdcCoefficients<_sdc_nodes,_multirate_nodes>::type
{
    
    enum
    {
        sdc_nodes = _sdc_nodes,
        multirate_nodes = _multirate_nodes,
        total_nodes = _total_nodes
    };
    value_type *Immk;
    SDCSpectralIntegrator<value_type,0,total_nodes,2,total_nodes> sdc_integrator;
    value_type coeff[sdc_nodes* ( total_nodes-1 ) ];
    value_type *points;
    value_type *dt;
    size_t ode_size;
    
    SDCSpectralIntegrator(size_t _ode_size) : Immk ( sdc_integrator.Immk ), points ( sdc_integrator.points ), dt ( sdc_integrator.dt ), ode_size(_ode_size)
    {
        get_coefficients ( coeff );
        value_type h = M_PI / ( sdc_nodes - 1 );
    }
    
    inline void integrate ( const value_type *F1, const value_type *F2 )
    {
        sdc_integrator.integrate ( F1 );
        for ( size_t j = 0; j < total_nodes - 1; ++j )
            for ( size_t i = 0; i < sdc_nodes; ++i )
                for ( size_t k = 0; k < ode_size; ++k )
                    Immk[j][k] += F2[i][k] * coeff[sdc_nodes*j+i];
    }
};


template<typename value_type,int _sdc_nodes>
struct SDCSpectralIntegrator<value_type, 0,_sdc_nodes,2,_sdc_nodes> : public SdcCoefficients<_sdc_nodes,0>::type
{
    
    enum
    {
        sdc_nodes = _sdc_nodes
    };
    
    value_type coeff[sdc_nodes* ( sdc_nodes-1 ) ];
    value_type points[sdc_nodes];
    value_type dt[sdc_nodes-1];
    value_type Dt;
    
    value_type* Immk[sdc_nodes-1];
    size_t ode_size;
    
    SDCSpectralIntegrator(size_t _ode_size) : ode_size(_ode_size)
    {
        for ( int i = 0; i < sdc_nodes-1; ++i )
            Immk[i] = new value_type[ode_size];
        get_coefficients ( coeff );
        value_type h = M_PI / ( sdc_nodes - 1 );
        for ( size_t i = 0; i < sdc_nodes; ++i )
        {
            value_type theta = -M_PI + i * h;
            points[i] = .5 * ( std::cos ( theta ) + 1 );
        }
        for ( size_t i = 0; i < sdc_nodes-1; ++i )
            dt[i] = points[i+1]-points[i];
    }
    ~SDCSpectralIntegrator()
    {
        for ( int i = 0; i < sdc_nodes-1; ++i )
            delete [] Immk[i];
    }
    inline void integrate ( value_type **F )
    {
        for ( size_t j = 0; j < sdc_nodes - 1; ++j )
        {
            std::fill ( Immk[j],Immk[j]+ode_size,value_type ( 0 ) );
            for ( size_t i = 0; i < sdc_nodes; ++i )
            {
                //                 std::copy(F[i],F[i]+ode_size,std::ostream_iterator<value_type>(std::cout, " "));
                for ( size_t k = 0; k < ode_size; ++k )
                    Immk[j][k] += F[i][k] * coeff[sdc_nodes*j+i];
            }
        }
    }
    
    inline void integrate ( value_type **F1, value_type **F2 )
    {
        value_type *F[sdc_nodes];
        for(size_t i = 0; i < sdc_nodes; ++i)
            F[i] = new value_type[ode_size];
        for ( size_t i = 0; i < sdc_nodes; ++i )
        {
            for ( size_t k = 0; k < ode_size; ++k )
                F[i][k] = F1[i][k]+F2[i][k];
            //             std::copy(F[i],F[i]+ode_size,std::ostream_iterator<value_type>(std::cout, " "));std::cout << std::endl;
        }
        integrate ( F );
        for(size_t i = 0; i < sdc_nodes; ++i)
            delete [] F[i];
    }
    inline void integrate ( const value_type **F1, const value_type **F2, const value_type **F3 )
    {
        value_type F[sdc_nodes][ode_size];
        for ( size_t i = 0; i < sdc_nodes; ++i )
            for ( size_t k = 0; k < ode_size; ++k )
                F[i][k] = F1[i][k]+F2[i][k];
            integrate ( F, F3 );
    }
    inline void integrate ( const value_type **F1, const value_type **F2, const value_type **F3, const value_type **F4 )
    {
        value_type F[sdc_nodes][ode_size];
        for ( size_t i = 0; i < sdc_nodes; ++i )
            for ( size_t k = 0; k < ode_size; ++k )
                F[i][k] = F1[i][k]+F2[i][k];
            integrate ( F, F3, F4 );
    }
    inline void integrate ( const value_type **F1, const value_type **F2, const value_type **F3, const value_type **F4, const value_type **F5 )
    {
        value_type F[sdc_nodes][ode_size];
        for ( size_t i = 0; i < sdc_nodes; ++i )
            for ( size_t k = 0; k < ode_size; ++k )
                F[i][k] = F1[i][k]+F2[i][k];
            integrate ( F, F3, F4, F5 );
    }
};

template <typename value_type, size_t _sdc_nodes = 3, size_t _ode_size = 0, size_t _multirate_nodes = 2, size_t _total_nodes = ( _sdc_nodes-1 ) * ( _multirate_nodes-1 ) + 1 >
struct ClenshawCurtis : public SDCSpectralIntegrator<value_type, _ode_size,_sdc_nodes,_multirate_nodes,_total_nodes>
    {};


#endif
