
#include<iostream>
#include<fstream>
#include<cstddef>
#include <cmath>
#include <iterator>
#include<algorithm>

#include "math/ode_solver/sdc/integrator/spectral_integrator.hpp"
#include "math/ode_solver/sdc/semi_implicit_sdc.hpp"
#include "rhs_functions.hpp"

int semi_implicit_sdc(int ac, char **av)
{
    typedef double value_type;


    {
	function_type<1> F;

        enum
        {
            sdc_nodes = 3,
            sdc_corrections = 3
        };
        const quadrature_type quad_name = gauss_lobatto;
        typedef Integrator<value_type, quad_name, sdc_nodes> spectral_integrator;
        SemiImplicitSDC<value_type, spectral_integrator, sdc_corrections > sdc(F.ode_size());

        const int size = 10;
        value_type x[size] = {1};
        value_type error[size] = {0};
        value_type time = 0, dt = .01;
        value_type f = 0;
        value_type x0 = x[0];
        int err = 0;
        for (size_t i = 0; i < size - 1; ++i)
        {
            x[i+1] = x[i];
            sdc(F, time, &x[i + 1], &f, dt);
            time += dt;
        }
    }

    {
        function_type<3> F;

        enum
        {
            sdc_nodes = 3,
            sdc_corrections = 3
        };
        const quadrature_type quad_type = gauss_lobatto;

        typedef Integrator<value_type, quad_type, sdc_nodes> spectral_integrator;
        SemiImplicitSDC<value_type, spectral_integrator, sdc_corrections > sdc(F.ode_size());

        const int size = 30;
        value_type x[size][3] = {{1, 0, 0}};
        value_type time = 0, dt = .01;
        value_type f[3];
        value_type x0[3] = {1, 0, 0};
        for (size_t i = 0; i < size - 1; ++i)
        {
            std::copy(x[i], x[i] + F.ode_size(), x[i+1]);
            sdc(F, time, x[i + 1], f, dt);
            time += dt;
        }
        x[0][0] = 1;x[0][1] = 0;x[0][2] = 0;
        std::cout << "x = [";
        for (int i = 0; i < size; ++i)
            std::copy(x[i], x[i] + 3, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "] \n";
    }
//
    {
        const int size = 100;
        const int spatial_size = 40;
        const int ode_size = 2 * spatial_size;

        diffusion_type<value_type, ode_size> F;
        enum
        {
            sdc_nodes = 3,
            sdc_corrections = 3
        };
        const quadrature_type quad_type = gauss_lobatto;

        typedef Integrator<value_type, quad_type, sdc_nodes> spectral_integrator;
        SemiImplicitSDC<value_type, spectral_integrator, sdc_corrections > sdc(F.ode_size());

        value_type **x = new value_type*[size];
        value_type **f = new value_type*[size];
        for (int i = 0; i < size; ++i)
        {
            x[i] = new value_type[ode_size];
            f[i] = new value_type[ode_size];
        }
        value_type time = 0, dt = .01;
        F.init(0, 1, x[0]);
        for (int i = 0; i < size - 1; ++i)
        {
            std::copy(x[i], x[i] + ode_size, x[i+1]);
            std::copy(f[i], f[i] + ode_size, f[i + 1]);
            sdc(F, time, x[i + 1], f[i + 1], dt);
            time += dt;
        }

        std::ofstream output("./data.m");
        output << "u = [";
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < ode_size; j += 2)
                output << x[i][j] << " ";
        output << "]; u = reshape(u," << spatial_size << "," << size << ")'; \n";
        output << "v = [";
        for (int i = 0; i < size; ++i)
            for (int j = 1; j < ode_size; j += 2)
                output << x[i][j] << " ";
        output << "]; v = reshape(v," << spatial_size << "," << size << ")';  \n";
        output << "[X Y] = meshgrid(linspace(0,1," << spatial_size << "),linspace(0," << size*dt << "," << size << "));\n";
        output << "mesh(X,Y,u)\n";
        output << "axis equal\n";
        for (int i = 0; i < size; ++i)
        {
            delete [] x[i];
            delete [] f[i];
        }
        delete [] x;
        delete [] f;
    }
}
