#ifndef SPRING_HPP
#define SPRING_HPP

#include<cmath>

/// Author: Ricardo Ortiz <ricardo.ortiz@tulane.edu>, (C) 2008
/// $Id: particle_system_spring.h 148 2010-02-18 20:53:38Z rortiz $

template<typename value_type, typename particle_type>
class Spring 
{

    private:
        particle_type * m_A;           ///< Pointer to one of the affected particles.
        particle_type * m_B;           ///< Pointer to one other affected particle.
        value_type      m_l;      ///< Rest length between the two particles.
        value_type      m_k;           ///< Spring Constant.


    public:

        particle_type       * A()              { return m_A; }
        particle_type       * B()              { return m_B; }
        particle_type const * A() const        { return m_A; }
        particle_type const * B() const        { return m_B; }
        value_type & stiffness()         { return m_k; }
        value_type const & stiffness()   const { return m_k; }
        value_type & resting_length()  { return m_l; }
        value_type const & resting_length() const { return m_l; }

    public:

        Spring() : m_A ( 0 ), m_B ( 0 ), m_l ( 0 ), m_k ( 1 ) { }
        ~Spring()  {  }

    public:

        /**
        * Init Spring.
        *
        * @param A
        * @param B
        */
        void init ( particle_type * A, particle_type * B )
        {
            assert ( A != B || !"Spring::init(): Particle A and B were the same" );
            assert ( A    || !"Spring::init(): Particle A was null" );
            assert ( B    || !"Spring::init(): Particle B was null" );

            this->m_A = A;
            this->m_B = B;
            
            value_type dx[3] = { A->position[0] - B->position[0], A->position[1] - B->position[1], A->position[2] - B->position[2] };
            m_l = std::sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
            
        }

    public:

        void apply()
        {
            value_type dx[3] = {m_A->position[0] - m_B->position[0],m_A->position[1] - m_B->position[1],m_A->position[2] - m_B->position[2]};
            value_type magnitude = std::sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
            assert ( magnitude > 0 || !"Spring::apply(): Non-positive spring length" );
            value_type L = m_k * ( 1.0 - m_l/magnitude );
            m_A->force[0] -= L*dx[0];
            m_A->force[1] -= L*dx[1];
            m_A->force[2] -= L*dx[2];
            m_B->force[0] += L*dx[0];
            m_B->force[1] += L*dx[1];
            m_B->force[2] += L*dx[2];
        }

};

#endif
