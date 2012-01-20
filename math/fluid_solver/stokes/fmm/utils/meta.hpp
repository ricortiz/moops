#ifndef FMM_META_HPP
#define FMM_META_HPP

#include "bits.hpp"

template<int i>
struct Factorial
{
    enum { value = i*Factorial<i-1>::value };
};

template<>
struct Factorial<1>
{
    enum { value = 1 };
};

template<>
struct Factorial<0>
{
    enum { value = 1 };
};

template<int n, int k>
struct Binomial
{
    enum
    {
        value = Factorial<n>::value/(Factorial<n-k>::value*Factorial<k>::value)
    };
};

template<size_t n>
struct Binomial<n,0>
{
    enum
    {
        value = 1
    };
};

template<size_t n>
struct Binomial<n,1>
{
    enum
    {
        value = n
    };
};

template<size_t n>
struct Binomial<n,n>
{
    enum
    {
        value = 1
    };
};


template<int dim>
struct inwhich_box 
{      
    template<typename T>
    static void result(const T *c, const T *p, int *tuple)
    {
        *tuple = (*p - *c < 0) ? 0 : 1;        
        inwhich_box<dim-1>::result(c+1,p+1,tuple+1);
    }
};

template<>
struct inwhich_box<0>
{      
    template<typename T>
    static void result(const T *c, const T *p, int *tuple)
    {
        *tuple = (*p - *c < 0) ? 0 : 1;
    }
};

#endif
