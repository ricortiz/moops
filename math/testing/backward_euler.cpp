
#include<iostream>
#include<cstddef>
#include <cmath>
#include <iterator>
#include<algorithm>

#include "math/ode_solver/euler/backward_euler.hpp"
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
    BackwardEuler<value_type,function_type<1> > backward_euler(F);
    
    const int size = 10000;
    
    value_type x[size] = {1}, rhs[size] = {0};
    value_type time = 0, dt = .001;
    for(size_t i = 0; i < size-1; ++i)
    {
        x[i+1] = x[i];
        backward_euler(time,&x[i+1],&rhs[i],dt);
        time += dt;
    }
    
    std::cout << "x = [";
    std::copy(x,x+size,std::ostream_iterator<value_type>(std::cout," "));
    std::cout << "]\n";
}