
#include<iostream>
#include<fstream>
#include<cstddef>
#include <cmath>
#include <iterator>
#include<algorithm>

#include "math/ode_solver/sdc/integrator/clenshaw_curtis.hpp"
#include "math/ode_solver/sdc/semi_implicit_sdc.hpp"
#include "rhs_functions.hpp"

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

    {
        const size_t size = 100;
        const size_t spatial_size = 40;
        const size_t ode_size = 2 * spatial_size;
        
        diffusion_type<value_type, ode_size> F;
        
        SemiImplicitSDC<value_type, diffusion_type<value_type,ode_size>, SDCSpectralIntegrator<value_type, 0, 5, 2, 5>, 5, 4 > sdc(F);
        
        value_type **x = new value_type*[size];
        value_type **f1 = new value_type*[size];
        value_type **f2 = new value_type*[size];
        for(size_t i = 0; i < size; ++i)
        {
            x[i] = new value_type[ode_size];
            f1[i] = new value_type[ode_size];
            f2[i] = new value_type[ode_size];
        }
        value_type time = 0, dt = .1;
        F.init(0, 1, x[0]);
        F.Explicit(time, x[0], f1[0]);
        F.Implicit(time, x[0], f2[0]);
        for(size_t i = 0; i < size - 1; ++i)
        {
            sdc(time, x[i + 1], x[i], f1[i + 1], f1[i], f2[i + 1], f2[i], dt);
            time += dt;
        }

        std::ofstream output("./data.m");
        output << "u = [";
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < ode_size; j += 2)
                output << x[i][j] << " ";
            output << "]; u = reshape(u," << spatial_size << "," << size << ")'; \n";
        output << "v = [";
        for(int i = 0; i < size; ++i)
            for(int j = 1; j < ode_size; j += 2)
                output<< x[i][j] << " ";
            output << "]; v = reshape(v," << spatial_size << "," << size << ")';  \n";
        output << "[X Y] = meshgrid(linspace(0,1," << spatial_size << "),linspace(0," << size*dt << "," << size << "));\n";
        output << "mesh(X,Y,u)\n";
        output << "axis equal\n";
        for(size_t i = 0; i < size; ++i)
        {
            delete [] x[i];
            delete [] f1[i];
            delete [] f2[i];
        }
        delete [] x;
        delete [] f1;
        delete [] f2;
    }
}
