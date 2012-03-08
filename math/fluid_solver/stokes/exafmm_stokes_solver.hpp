#ifndef EXAFMM_STOKES_SOLVER_HPP
#define EXAFMM_STOKES_SOLVER_HPP
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
#include <vector>
#include "external/exafmm/include/serialfmm.h"

template < typename value_type>
class ExaFmmStokesSolver
{
protected:
    typedef std::vector<Body>                   particle_array_type;
    typedef std::vector<Cell>                   box_array_type;
    typedef particle_array_type::iterator       particle_iterator;

    private:
        size_t                  m_num_particles;
        particle_array_type     m_particles;
        box_array_type          m_source_boxes;
        SerialFMM<Stokes>       m_fmm;
        bool m_images;

    public:
        ExaFmmStokesSolver(size_t num_particles)
                :
                m_num_particles(num_particles),
                m_particles(num_particles)                
        {
            m_fmm.initialize();            
        }

        box_array_type &getBoxes() { return m_source_boxes; }

        ~ExaFmmStokesSolver()
        {
            m_fmm.finalize();
        }
        
        inline void operator()(value_type, value_type *x, value_type *v, value_type *f)
        {
            size_t idx = 0;
            m_fmm.initTarget(m_particles);
            for( particle_iterator p = m_particles.begin(), end = m_particles.end(); p!= end; ++p )
            {
                p->X[0] = x[idx];
                p->X[1] = x[idx+1];
                p->X[2] = x[idx+2];
                
                p->FORCE[0] = f[idx]* 0.039788735772974;
                p->FORCE[1] = f[idx+1]* 0.039788735772974;
                p->FORCE[2] = f[idx+2]* 0.039788735772974;
                idx += 3;
            }
            m_fmm.setDomain(m_particles);
            m_source_boxes.clear();
            m_fmm.bottomup(m_particles, m_source_boxes);
            m_fmm.downward(m_source_boxes, m_source_boxes);
            for( particle_iterator p = m_particles.begin(), end = m_particles.end(); p!= end; ++p )
            {
                idx = 3*p->IBODY;
                v[idx] = p->TRG[0];
                v[idx+1] = p->TRG[1];
                v[idx+2] = p->TRG[2];
            }
//             idx = 0;
//             std::cout.precision(7);
//             std::cout << "velocities = [";
//             for( size_t i = 0; i < m_num_particles; ++i, idx+=3 )
//                 std::cout << v[idx] << " " << v[idx+1] << " " << v[idx+2] << " ";
//             std::cout << "];" << std::endl;
        }

        void setDelta(value_type delta) { m_fmm.setDelta(delta); }

        void allPairs()
        {
            m_fmm.buffer = m_particles;
            m_fmm.initTarget(m_fmm.buffer);
            m_fmm.evalP2P(m_fmm.buffer,m_fmm.buffer);
            real diff = 0, norm = 0;
            m_fmm.evalError(m_particles, m_fmm.buffer, diff, norm);
            m_fmm.printError(diff, norm);
        }
        SerialFMM<Stokes> &solver() { return m_fmm; }
        void withImages(bool images) { m_images = images; }
};


#endif


