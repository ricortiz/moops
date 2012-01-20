
#include<iostream>
#include<fstream>
#include<cstddef>
#include<cmath>
#include<iterator>

#include "math/ode_solver/euler/backward_euler.hpp"
#include "rhs_functions.hpp"

int main()
{
    typedef double value_type;
    {
        function_type<1> F;
        BackwardEuler<value_type, function_type<1> > backward_euler(F);

        const int size = 10;
        value_type error[size] = {0};

        value_type x[size] = {1}, rhs[size] = {0};
        value_type time = 0, dt = .001;
        for(size_t i = 0; i < size - 1; ++i)
        {
            x[i + 1] = x[i];
            backward_euler(time, &x[i + 1], &rhs[i], dt);
            time += dt;
	    error[i+1] = x[i+1] - std::exp(std::sin(time));
        }
        std::cout << "error = [";
        std::copy(error, error + size, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "]\n";
    }
//     {
//         function_type<3> F;
//         BackwardEuler<value_type, function_type<3> > backward_euler(F);
// 
//         const int size = 300;
// 
//         value_type x[size][3] = {{1, 0, 0}}, rhs[size][3];
//         value_type time = 0, dt = .001;
//         for(size_t i = 0; i < size - 1; ++i)
//         {
//             F.copy(x[i], x[i + 1]);
//             backward_euler(time, x[i + 1], rhs[i], dt);
//             time += dt;
//         }
//         std::cout << "x = [";
//         for(int i = 0; i < size; ++i)
//         {
//             std::copy(x[i], x[i] + 3, std::ostream_iterator<value_type>(std::cout, " "));
// //             std::cout << " ";
//         }
//         std::cout << "] \n";
//     }
//     {
//         const size_t size = 100;
//         const size_t spatial_size = 40;
//         const size_t ode_size = 2 * spatial_size;
// 
//         diffusion_type<value_type, ode_size> F;
// 
//         BackwardEuler<value_type, diffusion_type<value_type, ode_size> > euler(F);
// 
//         value_type **x = new value_type*[size];
//         value_type **f = new value_type*[size];
//         for(size_t i = 0; i < size; ++i)
//         {
//             x[i] = new value_type[ode_size];
//             f[i] = new value_type[ode_size];
//         }
//         value_type time = 0, dt = .1;
//         F.init(0, 1, x[0]);
//         F(time, x[0], f[0]);
//         for(size_t i = 0; i < size - 1; ++i)
//         {
//             euler(time, x[i + 1], x[i], f[i + 1], f[i], dt);
//             time += dt;
//         }
//         std::ofstream output("./data.m");
//         output << "u = [";
//         for(int i = 0; i < size; ++i)
//             for(int j = 0; j < ode_size; j += 2)
//                 output << x[i][j] << " ";
//             output << "]; u = reshape(u," << spatial_size << "," << size << ")'; \n";
//         output << "v = [";
//         for(int i = 0; i < size; ++i)
//             for(int j = 1; j < ode_size; j += 2)
//                 output<< x[i][j] << " ";
//             output << "]; v = reshape(v," << spatial_size << "," << size << ")';  \n";
//         
//         output << "[X Y] = meshgrid(linspace(0,1," << spatial_size << "),linspace(0," << size*dt << "," << size << "));\n";
//         output << "mesh(X,Y,u)\n";
//         output << "axis equal\n";
//         for(size_t i = 0; i < size; ++i)
//         {
//             delete [] x[i];
//             delete [] f[i];
//         }
//         delete [] x;
//         delete [] f;
//     }
}
