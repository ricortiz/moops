
#include<iostream>
#include<cstddef>
#include <cmath>
#include <iterator>
#include<algorithm>

#include "math/ode_solver/euler/forward_euler.hpp"
#include "rhs_functions.hpp"

int forward_euler(int ac, char **av)
{
    typedef double value_type;
    {
        function_type<1> F;
        ForwardEuler forward_euler(F.ode_size());

        const int size = 10;
        value_type error[size] = {0};

        value_type x[size] = {1}, rhs[size] = {0};
        value_type time = 0, dt = .0001;
        for (size_t i = 0; i < size - 1; ++i)
        {
            forward_euler(F,time,&x[i+1], &x[i], &rhs[i], dt);
            time += dt;
            error[i+1] = x[i+1] - std::exp(std::sin(time));
        }
	std::cout << "x = [";
        std::copy(x, x + size, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "]\n";
        std::cout << "error = [";
        std::copy(error, error + size, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "]\n";
    }
    {
        function_type<3> F;
        ForwardEuler forward_euler(F.ode_size());

        const int size = 300;

        value_type x[size][3] = {{1, 0, 0}}, rhs[size][3] = {{0}};
        value_type time = 0, dt = .001;
        for (size_t i = 0; i < size - 1; ++i)
        {
            forward_euler(F,time,x[i+1],x[i],rhs[i],dt);
            time += dt;
        }
        std::cout << "x = [";
        for (int i = 0; i < size; ++i)
            std::copy(x[i], x[i] + 3, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "] \n";
    }
}
