#ifndef SDC_NEWTON_KRYLOV_INTEGRATOR_HPP
#define SDC_NEWTON_KRYLOV_INTEGRATOR_HPP

#include "math/nonlinear_solver/inexact_newton.hpp"
#include "sdc_integrator.hpp"

template<typename boundary_type>
class SdcNewtonKrylovIntegrator
{
protected:
    typedef typename SDCIntegrator<boundary_type> explicit_sdc_type;

private:
    size_t    m_ode_size;
    explicit_sdc_type  m_sdc;
public:

    SdcNewtonKrylovIntegrator(size_t ode_size) : m_ode_size(ode_size), m_sdc(ode_size) {}

    template<typename value_type>
    void integrate(value_type timestep)
    {
        logger.startTimer("SdcNewtonKrylovIntegrator");
        value_type time = derived().time();
        m_sdc.predictor(time, timestep);
        m_sdc.corrector(time, timestep);
        logger.stopTimer("SdcNewtonKrylovIntegrator");
    }

    
    
};


#endif
