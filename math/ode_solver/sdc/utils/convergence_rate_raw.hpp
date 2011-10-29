#ifndef CONVERGENCE_RATE_HPP
#define CONVERGENCE_RATE_HPP

template<typename sdc_method, typename value_type>
void time_stepper(sdc_method &sdc,value_type dt,int iterations)
{
    value_type time = 0;
    for (int i = 0; i < iterations; ++i)
    {
        sdc.predictor(time,dt);
        sdc.corrector(time,dt);
        time += dt;
    }
}


template<typename sdc_type, typename value_type>
void test_convergence_rate(sdc_type &sdc, const value_type *x0, const value_type *f0, value_type final_time, const int *time_steps, int num_steps)
{    
    const int ode_size = sdc_traits<sdc_type>::ode_size;
    const value_type ref_time = final_time;
    value_type dt[num_steps];
    value_type ref_sol[ode_size];
    for(int i = 0; i < num_steps; ++i)
        dt[i] = final_time / time_steps[i];
        
    std::copy(x0,x0+ode_size,sdc.X(0));
    std::copy(f0,f0+ode_size,sdc.F(0));
    
    time_stepper(sdc,dt[0],time_steps[0]);
    std::copy(sdc.X(0),sdc.X(0)+ode_size,ref_sol);
    
    value_type error[num_steps];
    for(int i = 0; i < num_steps; ++i)
    {
        std::copy(x0,x0+ode_size,sdc.X(0));
        std::copy(f0,f0+ode_size,sdc.F(0));
        time_stepper(sdc,dt[i],time_steps[i]);
        error[i] = error_norm(ref_sol,sdc.X(0),ode_size);
    }
    std::copy(dt+1,dt+num_steps,std::ostream_iterator<value_type>(std::cout, " "));
    std::cout << std::endl;
    std::copy(error+1,error+num_steps,std::ostream_iterator<value_type>(std::cout, " "));
    std::cout << std::endl;

}

template<typename sdc_type, typename value_type>
void test_convergence_rate(sdc_type &sdc, const value_type *x0, const value_type *f1, const value_type *f2, value_type final_time, const int time_steps[], int num_steps)
{    
    const int ode_size = sdc_traits<sdc_type>::ode_size;
    const value_type ref_time = final_time;
    value_type dt[num_steps];
    value_type ref_sol[ode_size];
    for(int i = 0; i < num_steps; ++i)
        dt[i] = final_time / time_steps[i];
    
    sdc.X(0) = x0;
    sdc.Fi(0) = f1;
    sdc.Fe(0) = f2;
    time_stepper(sdc,dt[0],time_steps[0]);
    std::copy(sdc.X(0),sdc.X(0)+ode_size,ref_sol);
    
    value_type error[num_steps];
    for(int i = 1; i < num_steps; ++i)
    {        
        sdc.X(0) = x0;
        sdc.Fi(0) = f1;
        sdc.Fe(0) = f2;
        time_stepper(sdc,dt[i],time_steps[i]);
        error[i] = error_norm(ref_sol,sdc.X(0),ode_size);
    }
    std::copy(dt+1,dt+num_steps,std::ostream_iterator<value_type>(std::cout, " "));
    std::cout << std::endl;
    std::copy(error+1,error+num_steps,std::ostream_iterator<value_type>(std::cout, " "));
    std::cout << std::endl;
    
}

template<typename value_type>
value_type error_norm(value_type *a, value_type *b, size_t size)
{
    value_type s = value_type(0);
    for(size_t i = 0; i < size; ++i)
        s += (a[i]-b[i])*(a[i]-b[i]);
    return std::sqrt(s);    
}

#endif
