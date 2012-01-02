
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

    const size_t ode_size() {
        return 1;
    }
};

/// Robertson
template<>
struct function_type<3>
{
    template<typename value_type>
    void operator()(value_type t, const value_type *x, value_type *v)
    {
        v[0] = -.04*x[0] + 1e4*x[1]*x[2];
        v[1] = .04*x[0] - 1e4*x[1]*x[2]-3*1e7*x[1]*x[1];
        v[2] = 3*1e7*x[1]*x[1];
    }
    template<typename value_type>
    void copy(const value_type *x, value_type *y)
    {
        y[0] = x[0];
        y[1] = x[1];
        y[2] = x[2];
    }
    const size_t ode_size() {
        return 3;
    }
};

/// Diffusion (Brusselator)
template<typename value_type,size_t size>
struct difusion_type
{

    value_type dt;
    static const value_type alpha = 0.02;
    static const value_type B = 3.0;
    static const value_type A = 1.0;

    void operator()(value_type t, const value_type *y, value_type *v)
    {
        for (size_t i = 0; i < size; i+=2)
        {
            value_type y0 = (i == 0 ? 1.0 : y[i-2]);
            value_type yN = (i == size-2 ? 1.0 : y[i+2]);
            v[i] = A + y[i]*y[i]*y[i+1] - (B+1)*y[i] + alpha*1/dt*(y0-2*y[i]+yN);
        }
        for (size_t i = 1; i < size; i+=2)
        {
            value_type y0 = (i == 1 ? 3.0 : y[i-2]);
            value_type yN = (i == size-1 ? 3.0 : y[i+2]);
            v[i] = B*y[i-1] - y[i-1]*y[i-1]*y[i] + alpha*1.0/dt*(y0-2*y[i]+yN);
        }
    }
    void copy(const value_type *x, value_type *y)
    {
        std::copy(x,x+size,y);
    }

    void init(value_type x0, value_type xN, value_type *u0)
    {
        std::vector<value_type> x;
        size_t N = size/2;
        dt = (xN-x0)/N;
        x.resize(N);
        for (size_t i = 1; i < N+1; ++i)
            x[i-1] = value_type(i-1)/N;

        for (size_t i = 0, idx = 0; i < size; i+=2,++idx)
            u0[i] = 1.0 + std::sin(2*M_PI*x[idx]);
        for (size_t i = 1; i < size; i+=2)
            u0[i] = 3.0;

    }
    const size_t ode_size() {
        return size;
    }
};

int main()
{
    typedef double value_type;
    {
        function_type<1> F;

        ExplicitSDC<value_type,function_type<1>,SDCSpectralIntegrator<value_type, 0, 5, 2, 5>,5,4 > sdc(F);

        const int size = 10;

        value_type x[size] = {1}, rhs[size] = {0};
        value_type time = 0, dt = .1;
        value_type f0 = 1.0;
        value_type x0 = x[0];
        sdc.init(&x0,&f0);
        for (size_t i = 0; i < size-1; ++i)
        {
            sdc.predictor(time,dt);
            sdc.corrector(time,dt);
            time += dt;
            x[i+1] = sdc.X(0)[0]-std::exp(std::sin(time));
        }
        x[0] = 0;
        std::cout << "error = [";
        std::copy(x,x+size,std::ostream_iterator<value_type>(std::cout," "));
        std::cout << "]\n";
    }
    {
        function_type<3> F;

        ExplicitSDC<value_type, function_type<3>, SDCSpectralIntegrator<value_type, 0, 5, 2, 5>, 5, 4 > sdc(F);

        const int size = 30;
        value_type x[size][3] = {{1,0,0}};
        value_type time = 0, dt = .01;
        value_type f0[3];
        value_type x0[3] = {1, 0, 0};
        F(time,x0,f0);
        sdc.init(x0, f0);
        for (size_t i = 0; i < size - 1; ++i)
        {
            sdc.predictor(time, dt);
            sdc.corrector(time, dt);
            time += dt;
            F.copy(sdc.X(0),x[i+1]);
        }
        x[0][0] = 1;
        x[0][1] = 0;
        x[0][2] = 0;
        std::cout << "x = [";
        for (int i = 0; i < size; ++i)
            std::copy(x[i], x[i] + 3, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "] \n";
    }

    {
        const size_t size = 1000;
        const size_t spatial_size = 40;
        const size_t ode_size = 2*spatial_size;

        difusion_type<value_type,ode_size> F;

        ExplicitSDC<value_type, difusion_type<value_type,ode_size>, SDCSpectralIntegrator<value_type, 0, 5, 2, 5>, 5, 4 > sdc(F);

        value_type x[size][ode_size] = {0.0};
        value_type f[size][ode_size] = {0.0};
        value_type time = 0, dt = .01;
        F.init(0,1,x[0]);
        F(time,x[0],f[0]);
        sdc.init(x[0], f[0]);
        for (size_t i = 0; i < size - 1; ++i)
        {
            sdc.predictor(time, dt);
            sdc.corrector(time, dt);
            time += dt;
            F.copy(x[0],x[i+1]);
        }
        sdc.init(x[0], f[0]);
        std::cout << "u = [";
        for (int i = 0; i < size; ++i)
	  for(int j = 0; j < ode_size; j+=2)
            std::cout << x[i][j] << " ";
        std::cout << "] \n";
        std::cout << "v = [";
        for (int i = 0; i < size; ++i)
	  for(int j = 1; j < ode_size; j+=2)
            std::cout << x[i][j] << " ";
        std::cout << "] \n";
    }
}

