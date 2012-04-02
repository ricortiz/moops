#include<iostream>
#include<cassert>
#include<list>
#include<vector>
#include<map>
#include<fstream>

template<typename T>
struct Traits;

#include "particle_system/particle.hpp"
#include "particle_system/elastic_system/spring_system.hpp"

template<>
struct Traits<void>
{
    typedef ParticleWrapper<double> particle_type;
    typedef double value_type;
};


int spring_system(int, char **)
{
    SpringSystem<void> springs;
    Traits<void>::value_type x[3][3] = {{1,0,0},{0,1,0},{0,0,1}}, f[3][3] = {{0}};
    std::cout << "A = [" << x[0][0] << " " <<  x[0][1] << " " << x[0][2] << "];" << std::endl;
    std::cout << "B = [" << x[1][0] << " " <<  x[1][1] << " " << x[1][2] << "];" << std::endl;
    std::cout << "C = [" << x[2][0] << " " <<  x[2][1] << " " << x[2][2] << "];" << std::endl;
    Traits<void>::particle_type A, B, C;
    A.position = x[0]; A.force = f[0];
    B.position = x[1]; B.force = f[1];
    C.position = x[2]; C.force = f[2];

    Traits<void>::value_type k = 1.0;
    SpringSystem<void>::spring_iterator s;
    s = springs.addSpring(&A,&B,k); s->getAidx() = 0; s->getBidx() = 3;
    s = springs.addSpring(&A,&C,k); s->getAidx() = 0; s->getBidx() = 6;
    s = springs.addSpring(&B,&C,k); s->getAidx() = 3; s->getBidx() = 6;

    std::cout << springs;
    
    return 0;
}

