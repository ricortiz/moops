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
inline void computeStokeslet(const value_type *target, value_type *velocity, const value_type *source, const value_type *force, value_type delta)
{
    value_type dx[3] = {target[0] - source[0], target[1] - source[1], target[2] - source[2]};

    value_type r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    value_type d2 = delta * delta;
    value_type R1 = r2 + d2;
    value_type R2 = R1 + d2;
    value_type invR = 1.0 / R1;
    value_type H = std::sqrt(invR) * invR * 0.039788735772974;

    value_type fdx = (force[0] * dx[0] + force[1] * dx[1] + force[2] * dx[2]);

    velocity[0] += H * (force[0] * R2 + fdx * dx[0]);
    velocity[1] += H * (force[1] * R2 + fdx * dx[1]);
    velocity[2] += H * (force[2] * R2 + fdx * dx[2]);
}


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
inline void computeImage(const value_type *target, value_type *velocity, const value_type *source, const value_type *force, value_type delta)
{
    value_type h0 = source[2];

    value_type dx[3] = { target[0] - source[0],target[1] - source[1],target[2] + h0};
    
    value_type r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    value_type d2 = delta * delta;
    
    value_type R1 = r2 + d2;
    value_type R2 = R1 + d2;
    value_type invR = 1.0 / R1;
    value_type sqrtinvR = std::sqrt(invR);
    value_type sqrtinvR1 = sqrtinvR*invR;
    value_type sqrtinvR2 = sqrtinvR1*invR;

    value_type H1 = 0.5 * sqrtinvR*(1+d2*invR);
    value_type H2 = 0.5 * sqrtinvR1;

    value_type Hdp1 = -h0 * h0 * ( sqrtinvR1-3*d2*sqrtinvR2);
    value_type Hdp2 = -h0 * h0 * ( -3.0*sqrtinvR2 );

    value_type Hdb1 = h0 * ( -3.0 * d2 *sqrtinvR2 - sqrtinvR1 );
    value_type Hdb2 = h0 * ( sqrtinvR1 );
    value_type Hdb3 = h0 * ( -3.0 * sqrtinvR2 );

    value_type Hr1 = h0 * 3 * d2 * sqrtinvR2;

    // 21 flops
    value_type A11 = -(H1 + dx[0] * dx[0] * H2)*force[0];
    value_type A12 = -(dx[0] * dx[1] * H2)*force[1];
    value_type A13 = -(dx[0] * dx[2] * H2)*force[2];

    value_type A21 = -(dx[1] * dx[0] * H2)*force[0];
    value_type A22 = -(H1 + dx[1] * dx[1] * H2)*force[1];
    value_type A23 = -(dx[1] * dx[2] * H2)*force[2];

    value_type A31 = -(dx[2] * dx[0] * H2)*force[0];
    value_type A32 = -(dx[2] * dx[1] * H2)*force[1];
    value_type A33 = -(H1 + dx[2] * dx[2] * H2)*force[2];
    
    value_type fdx = (force[0] * dx[0] + force[1] * dx[1] + force[2] * dx[2]);
    
    value_type a = -H2 * (force[0] * R2 + fdx * dx[0]);
    value_type b = -H2 * (force[1] * R2 + fdx * dx[1]);
    value_type c = -H2 * (force[2] * R2 + fdx * dx[2]);

    //% A pos Dipole
    A11 = A11 - ( Hdp1 + dx[0] * dx[0] * Hdp2 )*force[0]; // %i=1 j=1
    A12 = A12 - ( dx[0] * dx[1] * Hdp2 )*force[1];    //  %i=1 j=2
    A13 = A13 + ( dx[0] * dx[2] * Hdp2 )*force[2];    //  %i=1 j=3

    A21 = A21 - ( dx[1] * dx[0] * Hdp2 )*force[0];    // %i=2 j=1
    A22 = A22 - ( Hdp1 + dx[1] * dx[1] * Hdp2 )*force[1]; //%i=2 j=2
    A23 = A23 + ( dx[1] * dx[2] * Hdp2 )*force[2];   //  %i=2 j=3

    A31 = A31 - ( dx[2] * dx[0] * Hdp2 )*force[0];   //  %i=3 j=1
    A32 = A32 - ( dx[2] * dx[1] * Hdp2 )*force[1];   //  %i=3 j=2
    A33 = A33 + ( Hdp1 + dx[2] * dx[2] * Hdp2 )*force[2]; //%i=3 j=3
    
    fdx = (-force[0] * dx[0] - force[1] * dx[1] + force[2] * dx[2]);
    a += (-force[0] * Hdp1+ fdx * dx[0] *Hdp2);
    b += (-force[1] * Hdp1+ fdx * dx[1] *Hdp2);
    c += (force[2] * Hdp1+ fdx * dx[2] *Hdp2);

    //% A pos Doublet
    A11 = A11 - dx[2] * (Hdb2 + dx[0] * dx[0] * Hdb3)*force[0] ;
    A12 = A12 - dx[2] *( dx[0] * dx[1] *  Hdb3 )*force[1];     //% i=1 j=2
    A13 = A13 + (dx[0] * Hdb2 + dx[0] * dx[2] *dx[2] * Hdb3 )*force[2];  // % i=1 j=3

    A21 = A21 - dx[2] * ( dx[1] * dx[0] *  Hdb3 )*force[0] ;  // % i=2 j=1
    A22 = A22 - dx[2] * (Hdb2 + dx[1] * dx[1] * Hdb3)*force[1] ;
    A23 = A23 + dx[1] * (Hdb2 + dx[2] * dx[2] * Hdb3 )*force[2] ; //% i=2 j=3

    A31 = A31 - dx[0] * (Hdb1 + dx[2] * dx[2] * Hdb3 )*force[0];
    A32 = A32 - dx[1] * (Hdb1 + dx[2] * dx[2] * Hdb3 )*force[1];
    A33 = A33 + dx[2] * (Hdb1 + 2 * Hdb2 + dx[2] * dx[2] * Hdb3)*force[2];

    fdx = (-force[0] * dx[0] - force[1] * dx[1] + force[2] * dx[2]);
    a += (-force[0] * Hdp1+ fdx * dx[0] *Hdp2);
    b += (-force[1] * Hdp1+ fdx * dx[1] *Hdp2);
    c += (force[2] * Hdp1+ fdx * dx[2] *Hdp2);
    
    //% A pos Rotlet
    A11 = A11 + dx[2] * Hr1*force[0];
    A22 = A22 + dx[2] * Hr1*force[1];
    A31 = A31 - dx[0] * Hr1*force[0];
    A32 = A32 - dx[1] * Hr1*force[1];
    velocity[0] += (A11 + A12 + A13  ) * 7.957747154594767e-02;
    velocity[1] += (A21  + A22  + A23 ) * 7.957747154594767e-02;
    velocity[2] += (A31  + A32  + A33 ) * 7.957747154594767e-02;
}

#endif
