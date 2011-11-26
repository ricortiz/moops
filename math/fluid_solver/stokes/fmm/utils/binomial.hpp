#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP

#include<sstream>
#include<cmath>
#include<stdexcept>
#include<typeinfo>

#include"factorial.hpp"

template <class T>
T binomial_coefficient(unsigned n, unsigned k)
{
    static const char* function = "binomial<%1%>(unsigned, unsigned)";
    if (k > n)
        return 0;
    T result;
    if ((k == 0) || (k == n))
        return 1;
    if ((k == 1) || (k == n-1))
        return n;

    if (n <= max_factorial<T>::value)
    {
        // Use fast table lookup:
        result = factorial<T>(n);
        result /= (factorial<T>(n-k)*factorial<T>(k));
    }
    else 
    {
        std::stringstream s;
        s << function << ": The binomial coefficient is undefined for n > " << max_factorial<T>::value;
        throw std::domain_error(s.str());
    }
    return result;
    // convert to nearest integer:
//     return std::ceil(result - 0.5f);
}


#endif
