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
        size_t m_num_sources;
        bool m_images;
        std::vector<value_type> buffer[2];
        
    public:
        idx_vector      col_ptr;
        idx_vector      col_idx;

    public:
        GpuStokesSolver(size_t num_sources) : m_num_sources(num_sources), m_images(false) {}
        
        inline void operator()(value_type, value_type *x, value_type *v, value_type *f)
        {
            operator()(0,x,v,x,f,m_num_sources);
        }

        inline void operator()(value_type, value_type *x, value_type *v, value_type *y, value_type *f, size_t num_targets)
        {
            std::fill(v,v+num_targets*3,0.0);
            ComputeStokeslets(x,v,y,f,m_delta,m_num_sources,num_targets,m_images);
        }

        inline void Implicit ( value_type, const value_type *x, value_type *v, const value_type *y, const value_type *f )
        {
            std::fill(v,v+m_num_sources*3,0.0);
            for ( size_t p = 0; p < col_ptr.size()-1; ++p )
                for(size_t j = col_ptr[p], end = col_ptr[p+1]; j != end; ++j)
                    ;
            
        }

        inline void Explicit ( value_type, const value_type *x, value_type *v, const value_type *y, const value_type *f )
        {
            std::fill(v,v+m_num_sources*3,0.0);
            std::vector<size_t>::iterator begin,end;
            for ( size_t p = 0; p < col_ptr.size()-1; ++p)
            {
                size_t sidx = 3*p;
                begin = col_idx.begin()+col_ptr[p];
                end = col_idx.begin()+col_ptr[p+1];
                // NOTE: col_idx sould have sorted column indices
                size_t min = col_idx[col_ptr[p]];
                size_t max = col_idx[col_ptr[p+1]-1];
                
                for(size_t j = 0; j < min; ++j)
                {
                    if (j == p) continue;
                    size_t tidx = 3*j;
                    computeStokeslet ( &x[tidx], &v[tidx], &y[sidx], &f[sidx], m_delta );
                }
                
                std::vector<size_t> tmp(max-min+1), idx(max-min+1-(col_ptr[p+1]-col_ptr[p]));
                generateSequence(tmp,min);
                std::set_difference(tmp.begin(),tmp.end(),begin,end,idx.begin());
                for(size_t j = 0; j < idx.size(); ++j)
                {
                    if (idx[j] == p) continue;
                    size_t tidx = 3*idx[j];
                    computeStokeslet ( &x[tidx], &v[tidx], &y[sidx], &f[sidx], m_delta );
                }
                
                for(size_t j = max+1; j < m_num_sources; ++j)
                {
                    if (j == p) continue;
                    size_t tidx = 3*j;
                    computeStokeslet ( &x[tidx], &v[tidx], &y[sidx], &f[sidx], m_delta );
                }
            }
        }
        
        void generateSequence(std::vector<size_t> &seq, size_t start)
        {
            for(size_t k = 0; k != seq.size(); ++k)
                seq[k] = start++;
        }
        
        void setDelta(value_type delta) { m_delta = delta; }
        void withImages(bool images) { m_images = images; }
        
};


#endif

