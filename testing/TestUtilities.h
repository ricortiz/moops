#ifndef TEST_UTILITIES_HPP
#define TEST_UTILITIES_HPP

#include<iostream>
#include<algorithm>
#include<limits>
#include <ctime>
#include <cmath>

template<typename value_type>
value_type dot(value_type *x,value_type *y, int size)
{
    value_type s = value_type(0);
    for (size_t i = 0; i < size; ++i)
        s += x[i]*y[i];
    return s;
}
template<typename value_type>
value_type error_norm(value_type *x,value_type *y, int size)
{
    value_type s = value_type(0);
    value_type z[size];
    std::transform(x,x+size,y,z,std::minus<value_type>());
    return std::sqrt(dot(z,z,size));
}

template<typename T>
bool compare(T *a, T *b, size_t size)
{
    return error_norm(a,b,size) < std::numeric_limits<T>::epsilon();
}

template<typename T>
T Random() { return T(rand()/T(RAND_MAX+1)); }


#endif
