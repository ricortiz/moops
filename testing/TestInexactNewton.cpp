
#include<iostream>
#include "math/nonlinear_solver/inexact_newton.hpp"

struct function_type
{
    template<typename value_type>
    void operator() ( const value_type *x, value_type *Fx )
    {
        Fx[0] = x[0]-x[1];
        Fx[1] = x[1]*x[1]-1;
    }
    
};

int TestInexactNewton(int ac, char *av[])
{
    typedef double value_type;
    const int size = 2;
    InexactNewtonMethod<value_type,size> newton;
    function_type F;
    value_type x[size] = {-2,-2};
    newton ( F,x,1e-9f,1e-9f );
    std::cout << x[0] << std::endl;
    std::cout << x[1] << std::endl;
    return 0;
}