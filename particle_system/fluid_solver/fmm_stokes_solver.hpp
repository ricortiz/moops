#ifndef FMM_STOKES_SOLVER
#define FMM_STOKES_SOLVER

#include "math/fluid_solver/stokes/fmm/octree/octree.hpp"
#include "math/fluid_solver/stokes/fmm/multipole_taylor.hpp"

template<typename particle_system_type, size_t max_particles, size_t precision>
class FMMStokesSolver
{
    protected:
        typedef typename particle_system_type::value_type value_type;
        typedef typename particle_system_type::particle_type particle_type;
        typedef Octree<value_type,max_particles,particle_type> tree_type;
        typedef MultipoleTaylor<value_type,precision> fmm_type;
        
    private:
        value_type m_delta;
        enum
        {
            num_particles = particle_system_type::num_particles
        };


    public:
        void operator()(value_type *x, value_type *v, value_type *f)
        {
            tree_type tree(x,v,f,num_particles);
            tree.init();
//             std::cout << tree << std::endl;
            fmm_type fmm(tree.boxes().size());
            for (int i = 0; i < 4; ++i)
            {
                fmm.compute_far_field_expansions(tree.root(), i);
                fmm.compute_local_expansions(tree.root(), i);
            }

            fmm.compute_velocity(tree.root(), m_delta);
//             std::cout << "v = [";
//             std::copy(v,v+3*num_particles,std::ostream_iterator<value_type>(std::cout," "));
//             std::cout << "];" <<  std::endl;
        }

        value_type const &delta() const { return m_delta; }
        value_type &delta() { return m_delta; }
};


#endif

