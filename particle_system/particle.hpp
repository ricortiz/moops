#ifndef PARTICLE_HPP
#define PARTICLE_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    Particle
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================

/// @name Particle - Stores the particle state.
/// @section Description Particle is the wrapper for the data arrays.  
/// 	It contains three pointers that gets set up when the particle system 
/// 	class is instantiated.
/// @section See Also 
/// 	ParticleSystem particle_system_storage

template<typename _value_type>
struct Particle
{
  typedef _value_type value_type;
    value_type *position;
    value_type *velocity;
    value_type *force;
    size_t i;
    size_t j;
    size_t k;
};


#endif
