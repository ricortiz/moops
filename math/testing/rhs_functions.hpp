#ifndef RHS_FUNCTIONS_HPP
#define RHS_FUNCTIONS_HPP

#include<vector>
#include<algorithm> // for std::copy()

///
/// Each functor defines the righ hand size of the ODE and overloads the following methods:
///
///     void operator()(value_type t, value_type *x, value_type *v)
///     void Explicit(value_type t, const value_type *x, value_type *v)
///     void Implicit(value_type t, const value_type *x, value_type *v)
///     size_t ode_size()
///
/// The operator()() method is used by the explicit methods and
/// Explicit() and Implicit() is used by the semi-implicit sdc method
///

template<int size>
struct function_type;

// dummy
template<>
struct function_type<1>
{
    template<typename value_type>
    void operator()(value_type t, const value_type *x, value_type *v)
    {
        v[0] = x[0] * std::cos(t);
    }

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
    size_t ode_size() {return 1u;}
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
    size_t ode_size() {return 3;}
};


/// Diffusion (Brusselator)
template<typename value_type, size_t size>
struct diffusion_type
{
    value_type dx;
    static const value_type alpha;
    static const value_type B;
    static const value_type A;

    void operator()(value_type t, const value_type *y, value_type *v)
    {
        size_t i = 0, j = 1;
        v[i] = A        + y[i] * y[i] * y[j] - (B + 1) * y[i] + alpha * 1.0 / (dx * dx) * (1.0 - 2 * y[i] + y[i + 2]);
        v[j] = B * y[i] - y[i] * y[i] * y[j]                  + alpha * 1.0 / (dx * dx) * (3.0 - 2 * y[j] + y[j + 2]);
        for(i = 2, j = 3; j < size - 2; i += 2, j += 2)
        {
            v[i] = A        + y[i] * y[i] * y[j] - (B + 1) * y[i] + alpha * 1.0 / (dx * dx) * (y[i - 2] - 2 * y[i] + y[i + 2]);
            v[j] = B * y[i] - y[i] * y[i] * y[j]                  + alpha * 1.0 / (dx * dx) * (y[j - 2] - 2 * y[j] + y[j + 2]);
        }
        i = size - 2, j = size - 1;
        v[i] = A               + y[i] * y[i] * y[j] - (B + 1) * y[i] + alpha * 1.0 / (dx * dx) * (y[i - 2] - 2 * y[i] + 1.0);
        v[j] = B * y[i]        - y[i] * y[i] * y[j]                         + alpha * 1.0 / (dx * dx) * (y[j - 2] - 2 * y[j] + 3.0);
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

    void Explicit(value_type t, const value_type *y, value_type *v)
    {
        for(size_t i = 0, j = 1; j < size; i += 2, j += 2)
        {
            v[i] = A        + y[i] * y[i] * y[j] - (B + 1) * y[i];
            v[j] = B * y[i] - y[i] * y[i] * y[j]                 ;
        }
    }

    void Implicit(value_type t, const value_type *y, value_type *v)
    {
        size_t i = 0, j = 1;
        v[i] = alpha * 1.0 / (dx * dx) * (1.0 - 2 * y[i] + y[i + 2]);
        v[j] = alpha * 1.0 / (dx * dx) * (3.0 - 2 * y[j] + y[j + 2]);
        for(i = 2, j = 3; j < size - 2; i += 2, j += 2)
        {
            v[i] = alpha * 1.0 / (dx * dx) * (y[i - 2] - 2 * y[i] + y[i + 2]);
            v[j] = alpha * 1.0 / (dx * dx) * (y[j - 2] - 2 * y[j] + y[j + 2]);
        }
        i = size - 2, j = size - 1;
        v[i] = alpha * 1.0 / (dx * dx) * (y[i - 2] - 2 * y[i] + 1.0);
        v[j] = alpha * 1.0 / (dx * dx) * (y[j - 2] - 2 * y[j] + 3.0);
    }

    size_t ode_size() { return size; }
};

template<typename value_type, size_t size>
const value_type diffusion_type<value_type,size>::alpha = 0.02;

template<typename value_type, size_t size>
const value_type diffusion_type<value_type,size>::B = 3;

template<typename value_type, size_t size>
const value_type diffusion_type<value_type,size>::A = 1;

#endif
