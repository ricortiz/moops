#ifndef FLUID_SOLVER_HPP
#define FLUID_SOLVER_HPP

template<typename Derived>
class FluidSolver
{
    protected:
        typedef typename surface_traits<Derived>::fluid_solver_type fluid_solver_type;

    private:
        fluid_solver_type m_fluid_solver;

    public:
        FluidSolver(size_t num_particles) : m_fluid_solver(num_particles) {}
        
        inline Derived &derived()
        {
            return *static_cast<Derived*>(this);
        }

        template<typename value_type>
        inline void operator()(value_type t, value_type *x, value_type *v)
        {
            derived().computeForces(t);
            m_fluid_solver(t,x,v,derived().forces());
        }

        template<typename value_type>
        inline void Explicit(value_type t, value_type *x, value_type *v)
        {
            derived().computeForces(t);
            m_fluid_solver.Explicit(t,x,v,derived().forces());
        }

        template<typename value_type>
        inline void Implicit(value_type t, value_type *x, value_type *v)
        {
            derived().computeForces(t);
            m_fluid_solver.Implicit(t,x,v,derived().forces());
        }
        fluid_solver_type &fluid_solver() { return m_fluid_solver; }
};



#endif
