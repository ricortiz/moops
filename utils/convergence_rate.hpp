#ifndef CONVERGENCE_RATE_HPP
#define CONVERGENCE_RATE_HPP

#include "sdc_base.hpp"
#include <Eigen/Core>
#include <Eigen/Array>

template<typename sdc_type, typename value_type, typename vector_type>
void test_convergence_rate(sdc_type &sdc, const vector_type &x0, const vector_type &f0, value_type final_time, const int time_steps[], int num_steps)
{    
    const value_type ref_time = final_time;
    value_type dt[num_steps];
    for(int i = 0; i < num_steps; ++i)
        dt[i] = final_time / time_steps[i];
    
    sdc.X(0) = x0;
    sdc.F(0) = f0;
    time_stepper(sdc,dt[0],time_steps[0]);
    vector_type ref_sol = sdc.X(0);
    
    value_type error[num_steps];
    for(int i = 0; i < num_steps; ++i)
    {
        sdc.X(0) = x0;
        sdc.F(0) = f0;
        time_stepper(sdc,dt[i],time_steps[i]);
        error[i] = (ref_sol - sdc.X(0)).norm();
    }
    std::copy(dt+1,dt+num_steps,std::ostream_iterator<value_type>(std::cout, " "));
    std::cout << std::endl;
    std::copy(error+1,error+num_steps,std::ostream_iterator<value_type>(std::cout, " "));
    std::cout << std::endl;

}

template<typename sdc_type, typename value_type, typename vector_type>
void test_convergence_rate(sdc_type &sdc, const vector_type &x0, const vector_type &f1, const vector_type &f2, value_type final_time, const int time_steps[], int num_steps)
{    
    const value_type ref_time = final_time;
    value_type dt[num_steps];
    for(int i = 0; i < num_steps; ++i)
        dt[i] = final_time / time_steps[i];
    
    sdc.X(0) = x0;
    sdc.Fi(0) = f1;
    sdc.Fe(0) = f2;
    time_stepper(sdc,dt[0],time_steps[0]);
    vector_type ref_sol = sdc.X(0);
    
    value_type error[num_steps];
    for(int i = 1; i < num_steps; ++i)
    {        
        sdc.X(0) = x0;
        sdc.Fi(0) = f1;
        sdc.Fe(0) = f2;
        time_stepper(sdc,dt[i],time_steps[i]);
        error[i] = (ref_sol - sdc.X(0)).norm();
    }
    std::copy(dt+1,dt+num_steps,std::ostream_iterator<value_type>(std::cout, " "));
    std::cout << std::endl;
    std::copy(error+1,error+num_steps,std::ostream_iterator<value_type>(std::cout, " "));
    std::cout << std::endl;
    
}
#endif
