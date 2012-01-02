
#include<iostream>
#include<cstddef>
#include <cmath>
#include <iterator>
#include<algorithm>

#include "math/ode_solver/euler/forward_euler.hpp"
#include <boost/type_traits/detail/is_mem_fun_pointer_impl.hpp>

template<int size>
struct function_type;

template<>
struct function_type<1>
{
    template<typename value_type>
    void operator()(value_type t, value_type *x, value_type *v)
    {
        v[0] = x[0] * std::cos(t);
    }

    const size_t ode_size() {return 1;}
};
/// Robertson
template<>
struct function_type<3>
{
    template<typename value_type>
    void operator()(value_type t, const value_type *x, value_type *v)
    {
        v[0] = -.04 * x[0] + 1e4 * x[1] * x[2];
        v[1] = .04 * x[0] - 1e4 * x[1] * x[2] - 3 * 1e7 * x[1] * x[1];
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
        ForwardEuler<value_type, function_type<1> > forward_euler(F);

        const int size = 10000;

        value_type x[size] = {1}, rhs[size];
        value_type time = 0, dt = .001;
        for (size_t i = 0; i < size - 1; ++i)
        {
            x[i+1] = x[i];
            forward_euler(time, &x[i+1], &rhs[i], dt);
            time += dt;
        }

        std::cout << "x = [";
        std::copy(x, x + size, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "]\n";
    }
    {
        function_type<3> F;
        ForwardEuler<value_type, function_type<3> > forward_euler(F);

        const int size = 300;

        value_type x[size][3] = {{1, 0, 0}}, rhs[size][3];
        value_type time = 0, dt = .001;
        for (size_t i = 0; i < size - 1; ++i)
        {
            F.copy(x[i], x[i+1]);
            forward_euler(time, x[i+1], rhs[i], dt);
            time += dt;
        }
        std::cout << "x = [";
        for (int i = 0; i < size; ++i)
            std::copy(x[i], x[i] + 3, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "] \n";
    }
}
