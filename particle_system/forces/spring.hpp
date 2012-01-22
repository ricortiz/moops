#ifndef SPRING_HPP
#define SPRING_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    Spring
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================
/// @name Spring - Enforces spring constraint between two particles
/// @section Description The Spring class is an abstraction of the spring force
///                      between any two points.
/// @section See also also SpringSystem

#include<cmath>

template<typename particle_type>
class Spring
{
    protected:
        typedef typename particle_type::value_type value_type;
    private:
        particle_type * m_A;           ///< Pointer to one of the affected particles.
        particle_type * m_B;           ///< Pointer to one other affected particle.
        value_type      m_l;      ///< Rest length between the two particles.
        value_type      m_k;           ///< Spring Constant.
        size_t          m_Aidx;           ///< Particle A index in particle array
        size_t          m_Bidx;           ///< Particle B index in particle array


    public:

        particle_type       * A()              { return m_A; }
        particle_type       * B()              { return m_B; }
        particle_type const * A() const        { return m_A; }
        particle_type const * B() const        { return m_B; }
        value_type & stiffness()         { return m_k; }
        value_type const & stiffness()   const { return m_k; }
        value_type & resting_length()  { return m_l; }
        value_type const & resting_length() const { return m_l; }
        size_t const &getAidx() const { return m_Aidx; }
        size_t const &getBidx() const { return m_Bidx; }
        size_t &getAidx() { return m_Aidx; }
        size_t &getBidx() { return m_Bidx; }

    public:

        Spring() : m_A(0), m_B(0), m_l(0), m_k(1) { }
        ~Spring()  {  }

    public:

        /**
        * Init Spring.
        *
        * @param A
        * @param B
        */
        void init(particle_type * A, particle_type * B)
        {
            assert(A != B || !"Spring::init(): Particle A and B were the same");
            assert(A    || !"Spring::init(): Particle A was null");
            assert(B    || !"Spring::init(): Particle B was null");

            this->m_A = A;
            this->m_B = B;

            value_type dx[3] = { A->position[0] - B->position[0], A->position[1] - B->position[1], A->position[2] - B->position[2] };
            m_l = std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

        }

    public:

        void apply()
        {
            apply(m_A->position, m_B->position, m_A->force, m_B->force);
        }

        inline void apply(const value_type *x1, const value_type *x2)
        {
            apply(x1, x2, m_A->force, m_B->force);
        }

        inline void apply(const value_type *x1, const value_type *x2, value_type *f1, value_type *f2)
        {
            value_type dx[3] = {x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2]};
            value_type l = std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
            assert(l > 0 || !"Spring::apply(): Non-positive spring length");
            value_type L = m_k * (1.0 - m_l / l);
            f1[0] -= L * dx[0];
            f1[1] -= L * dx[1];
            f1[2] -= L * dx[2];
            f2[0] += L * dx[0];
            f2[1] += L * dx[1];
            f2[2] += L * dx[2];
        }

};

#endif
