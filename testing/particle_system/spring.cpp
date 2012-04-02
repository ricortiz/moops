#include<iostream>
#include<cassert>
#include<list>
#include<vector>
#include<map>
#include<fstream>

#include "particle_system/particle.hpp"
#include "particle_system/forces/spring.hpp"

int spring(int, char **)
{
    typedef double value_type;
    typedef ParticleWrapper<value_type> particle_type;
    Spring<particle_type> spring;
    value_type x[2][3] = {{1,0,0},{0,1,0}}, f[2][3] = {{0}};
    std::cout << "A = [" << x[0][0] << " " <<  x[0][1] << " " << x[0][2] << "];" << std::endl;
    std::cout << "B = [" << x[1][0] << " " <<  x[1][1] << " " << x[1][2] << "];" << std::endl;
    particle_type A, B;
    A.position = x[0]; A.force = f[0];
    B.position = x[1]; B.force = f[1];

    spring.init(&A, &B);
    spring.getAidx() = 0;
    spring.getBidx() = 1;
    spring.stiffness() = 1;

    spring.apply();
    
    std::cout << spring;

    x[1][2] = 5;
    
    spring.apply();
    
    std::cout << spring;
    
    return 0;
}

