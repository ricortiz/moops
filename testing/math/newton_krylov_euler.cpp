
#include<iostream>
#include<cstddef>
#include <cmath>
#include <iterator>
#include<algorithm>

#include "math/ode_solver/euler/newton_krylov_euler.hpp"
#include "rhs_functions.hpp"

int main()
{
    typedef double value_type;
    {
        function_type<1> F;
        NewtonKrylovEuler<value_type, function_type<1> > euler(F);

        const int size = 10;
        value_type error[size] = {0};

        value_type x[size] = {1}, rhs[size];
        value_type time = 0, dt = .01;
        for (size_t i = 0; i < size - 1; ++i)
        {
            x[i+1] = x[i];
            euler(time, &x[i+1], &rhs[i], dt);
            time += dt;
            error[i+1] = x[i+1] - std::exp(std::sin(time));
        }

        std::cout << "error = [";
        std::copy(error, error + size, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "]\n";
    }
//     {
//         function_type<3> F;
//         ForwardEuler<value_type, function_type<3> > forward_euler(F);
// 
//         const int size = 300;
// 
//         value_type x[size][3] = {{1, 0, 0}}, rhs[size][3];
//         value_type time = 0, dt = .001;
//         for (size_t i = 0; i < size - 1; ++i)
//         {
//             F.copy(x[i], x[i+1]);
//             forward_euler(time, x[i+1], rhs[i], dt);
//             time += dt;
//         }
//         std::cout << "x = [";
//         for (int i = 0; i < size; ++i)
//             std::copy(x[i], x[i] + 3, std::ostream_iterator<value_type>(std::cout, " "));
//         std::cout << "] \n";
//     }
}
