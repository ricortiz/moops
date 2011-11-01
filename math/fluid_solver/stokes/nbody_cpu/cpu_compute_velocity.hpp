#ifndef CPU_COMPUTE_VELOCITY_HPP
#define CPU_COMPUTE_VELOCITY_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    compute_velocity
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
#include <cmath>

/**
 * @brief Update the velocity at target position due to a force at source position using a regularized stokeslet.
 *
 * @param target position of velocity
 * @param velocity velocity vector to update
 * @param source position of the force
 * @param force force vector
 * @param delta regularization parameter
 **/
template<typename value_type>
inline void compute_velocity(const value_type *target, value_type *velocity, const value_type *source, const value_type *force, value_type delta)
{
    value_type dx[3] = {target[0] - source[0],target[1] - source[1],target[2] - source[2]};
    
    value_type r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    value_type d2 = delta * delta;
    value_type R1 = r2 + d2;
    value_type R2 = R1 + d2;
    value_type invR = 1.0/R1;
    value_type H = std::sqrt(invR) *invR * 0.039788735772974;
    
    value_type fdx = (force[0]*dx[0]+force[1]*dx[1]+force[2]*dx[2]);
    
    velocity[0] += H* (force[0]*R2+fdx*dx[0]);
    velocity[1] += H* (force[1]*R2+fdx*dx[1]);
    velocity[2] += H* (force[2]*R2+fdx*dx[2]);
}

#endif
