#ifndef PARTICLE_HPP
#define PARTICLE_HPP

template<typename value_type>
struct Particle
{
    value_type *position;
    value_type *velocity;
    value_type *force;

// public:
//     Particle(): m_position(0), m_velocity(0), m_force(0) {}
//     inline value_type * position() { return m_position; }
//     inline value_type const * position() const { return m_position; }
//     inline value_type * velocity() { return m_velocity; }
//     inline value_type const * velocity() const { return m_velocity; }
//     inline value_type * force() { return m_force; }
//     inline value_type const * force() const { return m_force; }
//     friend void operator<< (std::ostream &out, const Particle<value_type> &p)
//     {
//         out << "--------------------------------" << std::endl;
//         out << "position = [" << p.m_position[0] << "," << p.m_position[1] << "," << p.m_position[2] << "]" << std::endl;
//         out << "velocity = [" << p.m_velocity[0] << "," << p.m_velocity[1] << "," << p.m_velocity[2] << "]" << std::endl;
//         out << "force = [" << p.m_force[0] << "," << p.m_force[1] << "," << p.m_force[2] << "]" << std::endl;
//         out << "--------------------------------" << std::endl;
//     }
};



#endif
