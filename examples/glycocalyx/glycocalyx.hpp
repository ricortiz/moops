#ifndef GLYCOCALYX_HPP
#define GLYCOCALYX_HPP

#include "particle_system/storage/particle_system_storage.hpp"
#include "particle_system/particle.hpp"
#include "particle_system/particle_system.hpp"
#include "geometry/tower_geometry.hpp"
#include "math/fluid_solver/stokes/cpu_stokes_solver.hpp"
#include "math/linear_solver/krylov/generalized_minimal_residual_method.hpp"

template<typename value_type, typename matrix_type>
struct MatrixOperator 
{
    const value_type *m_x;
    matrix_type &m_matrix;
    MatrixOperator (matrix_type &matrix): m_matrix(matrix){}
    void init(const value_type *x)
    {
        m_x = x;
    }
    
    inline void operator()(const value_type *y, value_type *Ay)
    {
        m_matrix(0,m_x,Ay,y);
    }
};

template<typename value_type, typename fluid_solver_type = CpuStokesSolver<value_type> >
class Glycocalyx : public ParticleSystem<Glycocalyx<value_type> >
{
    public:
        typedef ParticleSystem<Glycocalyx<value_type> > particle_system_type;

        typedef ParticleWrapper<value_type>             particle_type;
        typedef TowerGeometry<value_type>               tower_type;

    private:
        tower_type              m_geometry;             //< Contains geometrical props of sperms
        int                     m_num_geometries;
        fluid_solver_type       m_fluid_solver;
        MatrixOperator<value_type,fluid_solver_type> m_A;
        GeneralizedMinimalResidualMethod<value_type> m_linear_solver;

    public:

        Glycocalyx(size_t M, size_t N, int num_towers = 1)
        : particle_system_type(num_towers*M*N), m_fluid_solver(num_towers*M*N), m_num_geometries(num_towers), m_A(m_fluid_solver), m_linear_solver(3*num_towers*M*N)
        {
            // Set geometry parameters
            m_geometry.setDimensions(M, N);    // tail dims: MtxNt; head dims: MhxNh
            m_geometry.setLength(20.0);                   // length of the tail
            m_geometry.setRadius(1.0);                   // radius
            std::vector<value_type> mesh2d;              // coords of each geometry
            size_t x = std::sqrt(num_towers) , y = x;
            setGeometryGrid(mesh2d, x, y);         // create a grid to put the geometries on
            m_geometry.init(&this->particles() [0]);
            for(int i = 1, idx = 3; i < num_towers; ++i, idx += 3)
            {
                particle_type *p_init = &this->particles() [0];
                particle_type *p = &this->particles() [i * m_geometry.numParticles() ];
                for(size_t j = 0; j < m_geometry.numParticles(); ++j)
                {
                    p[j].position[0] = p_init[j].position[0] + mesh2d[idx];
                    p[j].position[1] = p_init[j].position[1] + mesh2d[idx + 1];
                    p[j].position[2] = p_init[j].position[2] + mesh2d[idx + 2];
                    p[j].i = p_init[j].i;
                    p[j].j = p_init[j].j;
                }
            }

        }

        inline void computeVelocities()
        {
            std::fill(this->velocities_begin(),this->velocities_end(),0.0);
            for(size_t i = 0, k = 2; i < this->data_size(); i+=3, k+=3)
                this->velocity()[i] = -.1*this->positions()[k];
            m_A.init(this->positions(),m_fluid_solver);
        }
        inline void computeForces()
        {
            m_linear_solver(m_A,this->velocities(),this->forces());
        }

        inline void run(value_type dt)
        {
            computeVelocities();
            computeForces();
        }


    private:
        void setGeometryGrid(std::vector<value_type> &mesh2d, int x, int y)
        {
            value_type dtheta = 1., dalpha = 1.;
            for(int i = 0; i < x; ++i)
                for(int j = 0; j < y; ++j)
                {
                    mesh2d.push_back(i * dtheta);
                    mesh2d.push_back(0.0);
                    mesh2d.push_back(j * dalpha);
                }
        }

        void getStrengths(const std::vector<size_t> &, const std::vector<size_t> &col_idx, std::vector<value_type> &strengths)
        {
            strengths.resize(col_idx.size(), 1.0);
        }

    public:
        tower_type &geometry()
        {
            return m_geometry;
        }
        
        fluid_solver_type &fluid_solver() { return m_flu}



};

template<typename _value_type>
struct Traits<Glycocalyx<_value_type> >
{
    typedef _value_type                                                 value_type;
    typedef ParticleWrapper<value_type>                                 particle_type;
    typedef ParticleSystemStorage<value_type, particle_type, SURFACE>   storage_type;
};

#endif

