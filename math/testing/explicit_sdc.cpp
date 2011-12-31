
#include<iostream>
#include<cstddef>
#include <cmath>
#include <iterator>
#include<algorithm>

#include "math/ode_solver/sdc/integrator/clenshaw_curtis.hpp"
#include "math/ode_solver/sdc/explicit_sdc.hpp"
#include <boost/type_traits/detail/is_mem_fun_pointer_impl.hpp>

template<int size>
struct function_type;

template<>
struct function_type<1>
{
    template<typename value_type>
    void operator()(value_type t, const value_type *x, value_type *v)
    {
        v[0] = x[0] * std::cos(t);
    }
    
    const size_t ode_size() {return 1;}
};

int main()
{
    typedef double value_type;
    function_type<1> F;
    
    ExplicitSDC<value_type,function_type<1>,SDCSpectralIntegrator<value_type, 0, 5, 2, 5>,5,4 > sdc(F);
    
    const int size = 10000;
    
    value_type x[size] = {1}, rhs[size] = {0};
    value_type time = 0, dt = .1;
    value_type f0 = 1.0;
    value_type x0 = x[0];
    sdc.init(&x[0],&f0);
    for(size_t i = 0; i < size-1; ++i)
    {
        sdc.predictor(time,dt);
        sdc.corrector(time,dt);
        x[i+1] = sdc.X(0)[0];
        time += dt;
    }
    
    std::cout << "x = [";
    std::copy(x,x+size,std::ostream_iterator<value_type>(std::cout," "));
    std::cout << "]\n";
}