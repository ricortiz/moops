#ifndef TOWER_GEOMETRY_H
#define TOWER_GEOMETRY_H

#include "geometry_base.h"

template<typename Types, int num_particles, int X, int Y>
class TowerGeometry : public Types::BaseGeometry<TowerGeometry<Types> >
{
public:
    typedef typename Types::value_type 	value_type;
    typedef typename Types::int_type 		int_type;
    typedef typename Types::vector3_type 	vector_3_type;


public:

    TowerGeometry(value_type length, value_type radius) {}




};

template<typename Types, int M, int N, int X, int Y>
struct geometry_traits<TowerGeometry<Types,_num_particles,N,M,X,Y> >
{
    typedef typename Types::value_type 		value_type;
    typedef typename Types::int_type 		int_type;
    typedef typename Types::vector3_type 	vector_3_type;
    enum
    {
        num_particles = M*N*X*Y,
        points_per_cross_section = M,
        num_cross_sections = N,
        grid_x = X,
        grid_y = Y
    };
};


#endif
