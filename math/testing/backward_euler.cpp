
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

/// Diffusion (Brusselator)
template<typename value_type, size_t size>
struct diffusion_type
{
    value_type dx;
    static const value_type alpha = 0.02;
    static const value_type B = 3.0;
    static const value_type A = 1.0;

    void operator()(value_type t, const value_type *y, value_type *v)
    {
        for(size_t i = 0; i < size; i += 2)
        {
            value_type y0 = (i == 0 ? 1.0 : y[i - 2]);
            value_type yN = (i == size - 2 ? 1.0 : y[i + 2]);
            v[i] = A + y[i] * y[i] * y[i + 1] - (B + 1) * y[i] + alpha * 1.0 / (dx * dx) * (y0 - 2 * y[i] + yN);
            assert(!isnan(v[i]));
        }
        for(size_t i = 1; i < size; i += 2)
        {
            value_type y0 = (i == 1 ? 3.0 : y[i - 2]);
            value_type yN = (i == size - 1 ? 3.0 : y[i + 2]);
            v[i] = B * y[i - 1] - y[i - 1] * y[i - 1] * y[i] + alpha * 1.0 / (dx * dx) * (y0 - 2 * y[i] + yN);
            assert(!isnan(v[i]));
        }
    }

    void copy(const value_type *x, value_type *y)
    {
        std::copy(x, x + size, y);
    }

    void init(value_type x0, value_type xN, value_type *v0)
    {
        std::vector<value_type> x;
        size_t N = size / 2;
        dx = (xN - x0) / N;
        x.resize(N);
        for(size_t i = 1; i < N + 1; ++i)
            x[i - 1] = value_type(i - 1) / N;

        for(size_t i = 0, idx = 0; i < size; i += 2, ++idx)
            v0[i] = 1.0 + std::sin(2 * M_PI * x[idx]);
        for(size_t i = 1; i < size; i += 2)
            v0[i] = 3.0;

    }
    const size_t ode_size()
    {
        return size;
    }
};

int main()
{
    typedef double value_type;
    {
        function_type<1> F;
        BackwardEuler<value_type, function_type<1> > backward_euler(F);

        const int size = 300;

        value_type x[size] = {1}, rhs[size] = {0};
        value_type time = 0, dt = .001;
        for(size_t i = 0; i < size - 1; ++i)
        {
            x[i + 1] = x[i];
            backward_euler(time, &x[i + 1], &rhs[i], dt);
            time += dt;
        }
        std::cout << "x = [";
        std::copy(x, x + size, std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "]\n";
    }
    {
        function_type<3> F;
        BackwardEuler<value_type, function_type<3> > backward_euler(F);

        const int size = 300;

        value_type x[size][3] = {{1, 0, 0}}, rhs[size][3];
        value_type time = 0, dt = .001;
        for(size_t i = 0; i < size - 1; ++i)
        {
            F.copy(x[i], x[i + 1]);
            backward_euler(time, x[i + 1], rhs[i], dt);
            time += dt;
        }
        std::cout << "x = [";
        for(int i = 0; i < size; ++i)
        {
            std::copy(x[i], x[i] + 3, std::ostream_iterator<value_type>(std::cout, " "));
//             std::cout << " ";
        }
        std::cout << "] \n";
    }
    {
        const size_t size = 1000;
        const size_t spatial_size = 40;
        const size_t ode_size = 2 * spatial_size;

        diffusion_type<value_type, ode_size> F;

        BackwardEuler<value_type, diffusion_type<value_type, ode_size> > euler(F);

        value_type **x = new value_type*[size];
        value_type **f = new value_type*[size];
        for(size_t i = 0; i < size; ++i)
        {
            x[i] = new value_type[ode_size];
            f[i] = new value_type[ode_size];
        }
        value_type time = 0, dt = .01*.5;
        F.init(0, 1, x[0]);
        F(time, x[0], f[0]);
        for(size_t i = 0; i < size - 1; ++i)
        {
            euler(time, x[i + 1], x[i], f[i + 1], f[i], dt);
            time += dt;
        }
        std::cout << "u = [";
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < ode_size; j += 2)
                std::cout << x[i][j] << " ";
        std::cout << "]; u = reshape(u," << spatial_size << "," << size << ")'; \n";
        std::cout << "v = [";
        for(int i = 0; i < size; ++i)
            for(int j = 1; j < ode_size; j += 2)
                std::cout << x[i][j] << " ";
        std::cout << "]; v = reshape(v," << spatial_size << "," << size << ")';  \n";
        for(size_t i = 0; i < size; ++i)
        {
            delete [] x[i];
            delete [] f[i];
        }
        delete [] x;
        delete [] f;
    }
}
