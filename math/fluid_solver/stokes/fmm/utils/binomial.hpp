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
//         throw std::domain_error(std::string(function) + std::string(": The binomial coefficient is undefined for k > n."));
    T result;
    if ((k == 0) || (k == n))
        return 1;
    if ((k == 1) || (k == n-1))
        return n;

    if (n <= max_factorial<T>::value)
    {
        // Use fast table lookup:
        result = factorial<T>(n);
        result /= factorial<T>(n-k);
        result /= factorial<T>(k);
    }
    else 
    {
        std::stringstream s;
        s << function << ": The binomial coefficient is undefined for n > " << max_factorial<T>::value;
        throw std::domain_error(s.str());
    }

    // convert to nearest integer:
    return std::ceil(result - 0.5f);
}


#endif
