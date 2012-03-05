#ifndef GPU_STOKES_SOLVER
#define GPU_STOKES_SOLVER
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

#include "math/fluid_solver/stokes/nbody_gpu/gpu_compute_velocity.hpp"

template<typename value_type>
class GpuStokesSolver
{        
    private:
        value_type m_delta;
        size_t m_num_particles;
        bool m_images;
        
    public:
        GpuStokesSolver(size_t num_particles) : m_num_particles(num_particles), m_images(false) {}
        
        inline void operator()(value_type, value_type *x, value_type *v, value_type *f)
        {
            operator()(0,x,v,x,f,m_num_particles,m_num_particles);
        }

        inline void operator()(value_type, value_type *x, value_type *v, value_type *y, value_type *f, size_t num_sources, size_t num_targets)
        {
            std::fill(v,v+num_targets*3,0.0);
            ComputeStokeslets(x,v,y,f,m_delta,num_sources,num_targets,m_images);
        }
        
        inline void operator()(value_type, value_type *x, value_type *v, value_type *y, value_type *f, size_t num_targets)
        {
            std::fill(v,v+num_targets*3,0.0);
            ComputeStokeslets(x,v,y,f,m_delta,m_num_particles,num_targets,m_images);
        }

        void setDelta(value_type delta) { m_delta = delta; }
        void withImages(bool images) { m_images = images; }
        
};


#endif

