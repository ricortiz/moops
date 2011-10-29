#ifndef EULER_INTEGRATOR_HPP
#define EULER_INTEGRATOR_HPP

template<typename T> struct immersed_structure_traits;

template<typename boundary_type>
class EulerIntegrator
{
    protected:
        typedef typename immersed_structure_traits<boundary_type>::value_type          value_type;
        typedef typename immersed_structure_traits<boundary_type>::particle_integrator_type integrator_type;

    private:
        integrator_type euler_integrator;

    public:

        inline boundary_type &boundary()
        {
            return *static_cast<boundary_type*>(this);
        }

        template<typename value_type>
        inline void integrate(value_type timestep)
        {
            value_type time = boundary().time();
            euler_integrator(boundary(),time,timestep);
        }

};



#endif
