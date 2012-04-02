#include<cassert>
#include<list>
#include<vector>
#include<map>
#include<fstream>
#include "utils/logger.hpp"

// Logger logger;

template<typename T>
struct Traits;

#include "particle_system/particle.hpp"
#include "particle_system/elastic_system/elastic_boundary.hpp"

template<>
struct Traits<int>
{
    typedef ParticleWrapper<double> particle_type;
    typedef double value_type;
};


int elastic_system(int, char **)
{
    ElasticBoundary<int> elastic_boundary;
    Traits<int>::value_type x[3][3] = {{1,0,0},{0,1,0},{0,0,1}}, f[3][3] = {{0}};
    std::cout << "A = [" << x[0][0] << " " <<  x[0][1] << " " << x[0][2] << "];" << std::endl;
    std::cout << "B = [" << x[1][0] << " " <<  x[1][1] << " " << x[1][2] << "];" << std::endl;
    std::cout << "C = [" << x[2][0] << " " <<  x[2][1] << " " << x[2][2] << "];" << std::endl;
    Traits<int>::particle_type A, B, C;
    A.position = x[0]; A.force = f[0];
    B.position = x[1]; B.force = f[1];
    C.position = x[2]; C.force = f[2];

    Traits<int>::value_type k = 1.0;
    ElasticBoundary<int>::spring_iterator s;
    s = elastic_boundary.addSpring(&A,&B,k); s->getAidx() = 0; s->getBidx() = 3;
    s = elastic_boundary.addSpring(&A,&C,k); s->getAidx() = 0; s->getBidx() = 6;
    s = elastic_boundary.addSpring(&B,&C,k); s->getAidx() = 3; s->getBidx() = 6;

    x[0][0] = 2; x[0][1] = 0; x[0][2] = 0;
    x[1][0] = 0; x[1][1] = 2; x[1][2] = 0;
    
    elastic_boundary.computeForces();

    std::cout << "fA = [" << f[0][0] << " " <<  f[0][1] << " " << f[0][2] << "];" << std::endl;
    std::cout << "fB = [" << f[1][0] << " " <<  f[1][1] << " " << f[1][2] << "];" << std::endl;
    std::cout << "fC = [" << f[2][0] << " " <<  f[2][1] << " " << f[2][2] << "];" << std::endl;
    
    x[0][0] = 1; x[0][1] = 0; x[0][2] = 0;
    x[1][0] = 0; x[1][1] = 1; x[1][2] = 0;

    std::fill(&f[0][0], &f[0][0]+9,0.0);
    Traits<int>::value_type y[9] = {2,0,0,0,2,0,0,0,1};
    elastic_boundary.computeForces(y);
    
    std::cout << "fA = [" << f[0][0] << " " <<  f[0][1] << " " << f[0][2] << "];" << std::endl;
    std::cout << "fB = [" << f[1][0] << " " <<  f[1][1] << " " << f[1][2] << "];" << std::endl;
    std::cout << "fC = [" << f[2][0] << " " <<  f[2][1] << " " << f[2][2] << "];" << std::endl;
    
    std::fill(&f[0][0], &f[0][0]+9,0.0);
    elastic_boundary.computeForces(y,&f[0][0]);
    
    std::cout << "fA = [" << f[0][0] << " " <<  f[0][1] << " " << f[0][2] << "];" << std::endl;
    std::cout << "fB = [" << f[1][0] << " " <<  f[1][1] << " " << f[1][2] << "];" << std::endl;
    std::cout << "fC = [" << f[2][0] << " " <<  f[2][1] << " " << f[2][2] << "];" << std::endl;

    
    
    return 0;
}

