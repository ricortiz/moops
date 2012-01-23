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
        inline Derived &derived()
        {
            return *static_cast<Derived*>(this);
        }

        template<typename value_type>
        inline void operator()(value_type t, value_type *x, value_type *v)
        {
            logger.startTimer("TimeStep");
            derived().updateForces(t);
            m_fluid_solver(t,x,v,derived().forces(),derived().particles_size());
            logger.stopTimer("TimeStep");
        }

        fluid_solver_type &fluid_solver() { return m_fluid_solver; }
};



#endif
