#ifndef FMM_STOKES_SOLVER
#define FMM_STOKES_SOLVER
/****************************************************************************
** MOOPS -- Modular Object Oriented Particle Simulator
** Copyright (C) 2011-2012  Ricardo Ortiz <ortiz@unc.edu>
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/
/// @name FMMStokesSolver - Wrapper for the Fast Multipole Method based Stokes solver
/// @section Description FMMStokesSolver wraps the FMM method so the simulator can use it.
/// @section See also MultipoleTaylor

#include "math/fluid_solver/stokes/fmm/octree/octree.hpp"
#include "math/fluid_solver/stokes/fmm/multipole_taylor.hpp"
#include "particle_system/particle.hpp"

template<typename value_type, size_t max_particles, size_t precision>
class FmmStokesSolver
{
    protected:
        typedef ParticleWrapper<value_type> particle_type;
        typedef Octree<value_type,max_particles,particle_type> tree_type;
        typedef MultipoleTaylor<value_type,precision> fmm_type;
        
    private:
        value_type m_delta;
        size_t m_num_particles;
        bool m_images;
        
    public:
        FmmStokesSolver(size_t num_particles) : m_num_particles(num_particles) {}
        
        inline void operator()(value_type, value_type *x, value_type *v, value_type *, value_type *f)
        {
            tree_type tree(x,v,f,m_num_particles);
            tree.init();
            fmm_type fmm(tree.boxes().size());
            for (int i = 0; i < 4; ++i)
            {
                fmm.compute_far_field_expansions(tree.root(), i);
                fmm.compute_local_expansions(tree.root(), i);
            }

            fmm.compute_velocity(tree.root(), m_delta);
        }
        
        void setDelta(value_type delta) { m_delta = delta; }
        void withImages(bool images) { m_images = images; }
};


#endif

