#ifndef CPU_STOKES_SOLVER_HPP
#define CPU_STOKES_SOLVER_HPP

#include <vector>
#include <algorithm>
#include<iterator>
#include "nbody_cpu/cpu_compute_velocity.hpp"

template<typename value_type>
class CpuStokesSolver
{
    protected:
        typedef std::vector<size_t> idx_vector;

    private:
        value_type      m_delta;
        size_t          m_num_particles;

    public:    
        idx_vector      col_ptr;
        idx_vector      col_idx;

    public:
        CpuStokesSolver(size_t num_particles) : m_num_particles(num_particles) {}
        
        inline void operator() ( value_type t, const value_type *x, value_type *v, const  value_type *f)
        {
            operator() ( t, x, v, x, f, m_num_particles, m_num_particles );
        }

        inline void operator() ( value_type, const value_type *x, value_type *v, const value_type *y, const value_type *f, size_t num_sources, size_t num_targets )
        {
            size_t size_targets = 3 * num_targets;
            size_t size_sources = 3 * num_sources;
            std::fill(v,v+size_targets,0.0);
            #pragma omp parallel for shared(x,v,y,f)
            for ( size_t i = 0; i < size_targets; i += 3 )
                for ( size_t j = 0; j < size_sources; j += 3 )
                    compute_velocity ( &x[i], &v[i], &y[j], &f[j], m_delta );
//             std::cout << "direct_velocities = [";std::copy(v, v + 3*m_num_particles, std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "]" << std::endl;
        }

        inline void Implicit( value_type t, const value_type *x, value_type *v, const value_type *f )
        {
            Implicit( t, x, v, x, f );
        }

        inline void Explicit(value_type t, const value_type *x, value_type *v, const value_type *f )
        {
            Explicit  ( t, x, v, x, f );
        }

        inline void Implicit ( value_type, const value_type *x, value_type *v, const value_type *y, const value_type *f )
        {
            std::fill(v,v+m_num_particles*3,0.0);
            #pragma omp parallel for shared(x,v,y,f)
            for ( size_t p = 0; p < col_ptr.size()-1; ++p )
            {
                size_t sidx = 3*p;
                compute_velocity ( &x[sidx], &v[sidx], &y[sidx], &f[sidx], m_delta );
                for(size_t j = col_ptr[p], end = col_ptr[p+1]; j != end; ++j)
                {
                    size_t tidx = 3*col_idx[j];
                    compute_velocity ( &x[tidx], &v[tidx], &y[sidx], &f[sidx], m_delta );
                }
            }
        }

        inline void Explicit ( value_type, const value_type *x, value_type *v, const value_type *y, const value_type *f )
        {
            std::fill(v,v+m_num_particles*3,0.0);
            std::vector<size_t>::iterator begin,end;
            #pragma omp parallel for shared(x,v,y,f) private(begin,end)
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
                    compute_velocity ( &x[tidx], &v[tidx], &y[sidx], &f[sidx], m_delta );
                }
                
                std::vector<size_t> tmp(max-min+1), idx(max-min+1-(col_ptr[p+1]-col_ptr[p]));
                generate_sequence(tmp,min);
                std::set_difference(tmp.begin(),tmp.end(),begin,end,idx.begin());
                for(size_t j = 0; j < idx.size(); ++j)
                {
                    if (idx[j] == p) continue;
                    size_t tidx = 3*idx[j];
                    compute_velocity ( &x[tidx], &v[tidx], &y[sidx], &f[sidx], m_delta );
                }
                
                for(size_t j = max+1; j < m_num_particles; ++j)
                {
                    if (j == p) continue;
                    size_t tidx = 3*j;
                    compute_velocity ( &x[tidx], &v[tidx], &y[sidx], &f[sidx], m_delta );
                }
            }
        }

        
        void generate_sequence(std::vector<size_t> &seq, size_t start)
        {
            for(size_t k = 0; k != seq.size(); ++k)
                seq[k] = start++;
        }
      
        void setDelta(value_type delta) { m_delta = delta; }
};


#endif

