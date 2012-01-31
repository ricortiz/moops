#include<iostream>
#include<limits>
#include<cmath>
#include<map>
#include "math/nonlinear_solver/inexact_newton.hpp"

template<int size>
struct function_type;

template<>
struct function_type<2>
{
    template<typename value_type>
    void operator() ( const value_type *x, value_type *Fx )
    {
        Fx[0] = x[0]-x[1]+std::cos(x[1])*std::sin(x[1])-std::exp(-x[1]);
        Fx[1] = -x[1]*x[1]+std::sin(x[0])*std::cos(x[0])*std::exp(-x[0]);
    }
    
};

template<typename stats_type>
void display_stats(stats_type &stats)
{
    typedef typename stats_type::mapped_type::value_type value_type;
    typedef typename stats_type::iterator stats_iterator;
    typedef typename std::vector<value_type>::iterator data_iterator;
    typedef typename stats_type::key_type key_type;
    
    for(stats_iterator s = stats.begin(); s != stats.end(); ++s)
    {
        key_type key = s->first;
        std::vector<value_type> data = s->second;
        std::cout << key << " = [";
        std::copy(data.begin(),data.end(),std::ostream_iterator<value_type>(std::cout, " "));
        std::cout << "];\n";
    }
}

int inexact_newton(int , char **)
{
    typedef double value_type;
    typedef std::string key_type;
    typedef std::map<key_type,std::vector<value_type> > stats_type;

    stats_type stats;
    const int size = 2;
    InexactNewtonMethod<value_type,1000,100> newton(size);
    function_type<2> F;
    value_type x[size] = {1,1};
    newton ( F,x,1e-16f,1e-13f, &stats );

    display_stats(stats);
    return 0;
    
}
