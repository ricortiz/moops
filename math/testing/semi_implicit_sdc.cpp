
#include<iostream>
#include<cstddef>
#include <cmath>
#include <iterator>
#include<algorithm>

#include "math/ode_solver/sdc/integrator/clenshaw_curtis.hpp"
#include "math/ode_solver/sdc/semi_implicit_sdc.hpp"
#include <boost/type_traits/detail/is_mem_fun_pointer_impl.hpp>

template<int size>
struct function_type;

template<>
struct function_type<1>
{
    template<typename value_type>
    void Explicit(value_type t, const value_type *x, value_type *v)
    {
        v[0] = .99999 * x[0] * std::cos(t);
    }
    template<typename value_type>
    void Implicit(value_type t, const value_type *x, value_type *v)
    {
        v[0] = .00001 * x[0] * std::cos(t);
    }
    const size_t ode_size() {return 1;}
};

/// Robertson
template<>
struct function_type<3>
{
    template<typename value_type>
    void Explicit(value_type t, const value_type *x, value_type *v)
    {
        v[0] = -.04 * x[0];
        v[1] = .04 * x[0];
        v[2] = 0.0;
    }
    template<typename value_type>
    void Implicit(value_type t, const value_type *x, value_type *v)
    {
        v[0] = 1e4 * x[1] * x[2];
        v[1] = - 1e4 * x[1] * x[2] - 3 * 1e7 * x[1] * x[1];
        v[2] = 3 * 1e7 * x[1] * x[1];
    }
    
    template<typename value_type>
    void copy(const value_type *x, value_type *y)
    {
        y[0] = x[0];
        y[1] = x[1];
        y[2] = x[2];
    }
    const size_t ode_size() {return 3;}
};

int main()
{
    typedef double value_type;


    {
        function_type<1> F;

        SemiImplicitSDC<value_type, function_type<1>, SDCSpectralIntegrator<value_type, 0, 5, 2, 5>, 5, 4 > sdc(F);

        const int size = 10000;
        value_type x[size] = {1};
        value_type time = 0, dt = .1;
        value_type f01 = .5, f02 = .5;
        value_type x0 = x[0];
        sdc.init(&x0, &f01, &f02);
        for (size_t i = 0; i < size - 1; ++i)
        {
            sdc.predictor(time, dt);
            sdc.corrector(time, dt);
            time += dt;
            x[i+1] = sdc.X(0)[0] - std::exp(std::sin(time));
        }
        x[0] = 0;
        std::cout << "error = [";
        std::copy(x, x + size, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "]\n";
    }

    {
        function_type<3> F;

        SemiImplicitSDC<value_type, function_type<3>, SDCSpectralIntegrator<value_type, 0, 5, 2, 5>, 5, 4 > sdc(F);

        const int size = 30;
        value_type x[size][3] = {{1,0,0}};
        value_type time = 0, dt = .01;
        value_type f01[3], f02[3];
        value_type x0[3] = {1, 0, 0};
	F.Implicit(time,x0,f01);
	F.Explicit(time,x0,f02);
        sdc.init(x0, f01, f02);
        for (size_t i = 0; i < size - 1; ++i)
        {
            sdc.predictor(time, dt);
            sdc.corrector(time, dt);
            time += dt;
            F.copy(sdc.X(0),x[i+1]);
        }
        x[0][0] = 1;x[0][1] = 0;x[0][2] = 0;
        std::cout << "x = [";
        for (int i = 0; i < size; ++i)
            std::copy(x[i], x[i] + 3, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "] \n";
    }
}
