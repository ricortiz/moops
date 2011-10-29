#ifndef SWIMMERS_HPP
#define SWIMMERS_HPP

#include <cassert>
#include "surface.hpp"
#include "sine_geometry.hpp"
#include <deque>


template<typename value_type>
class Swimmers : public Surface<Swimmers<value_type> >
{
    public:
        typedef SineGeometry<value_type> sine_type;
        typedef typename Surface<Swimmers<value_type> >::grid_type grid_type;

    private:
        BaseGeometry<sine_type> *m_swimmers;
        std::vector<size_t>    m_col_ptr;
        std::vector<size_t>    m_col_idx;

    public:

        void init(size_t num_swimmers, BaseGeometry<sine_type> &sine_geometry, value_type *positions, size_t num_particles)
        {
            typedef std::vector<size_t>::iterator iterator;
            size_t num_points = num_particles/num_swimmers;
            assert(num_points == sine_geometry.size());
            m_col_ptr.push_back(0);
            grid_type &grid = this->grid();
            size_t nx = std::sqrt(num_swimmers)+1;
            size_t ny = nx;
            std::deque<value_type> X;
            value_type dx[] = {1.5,.5};
            get_squared_patch(X,nx,ny,dx);
            for (size_t i = 0; i < num_swimmers; ++i)
            {
                value_type *x0 = sine_geometry.get_x0();
                x0[0] = X.front();
                X.pop_front();
                x0[1] = X.front();
                X.pop_front();
                size_t size = m_col_idx.size();
                size_t offset = i*num_points;
                sine_geometry.init(&positions[3*i*num_points],grid);
                sine_geometry.get_connections(m_col_ptr,m_col_idx);
                for (iterator i = m_col_idx.begin()+size, end = m_col_idx.end(); i != end; ++i)
                    *i += offset;
            }
//             for (size_t p = 0; p < m_col_ptr.size()-1; ++p)
//             {
//                 std::cout << "p[" << p << "] = [";
//                 for (size_t i = m_col_ptr[p], end = m_col_ptr[p+1]; i < end; ++i)
//                     std::cout << m_col_idx[i] << ",";
//                 std::cout << "]" << std::endl;
//             }
            m_swimmers = &sine_geometry;
        }

        template<typename spring_system_type>
        void set_springs(spring_system_type &spring_system)
        {
            typedef typename spring_system_type::particle_type particle_type;
            particle_type *particles = spring_system.particles();
            for (size_t p = 0; p < m_col_ptr.size()-1; ++p)
                for (size_t i = m_col_ptr[p], end = m_col_ptr[p+1]; i < end; ++i)
                    if (!spring_system.exist_spring(&particles[p],&particles[m_col_idx[i]]))
                        spring_system.add_spring(&particles[p],&particles[m_col_idx[i]]);
            std::cout << "Created " << spring_system.springs_size() << " springs." << std::endl;
//             int i = 0;
//             grid_type &grid = this->grid();
//             typedef typename spring_system_type::spring_iterator iterator;
//             for (iterator s = spring_system.spring_begin(), end = spring_system.spring_end(); s != end; ++s, ++i)
//                 std::cout << "spring[" << i << "] = [(" << grid[s->A()->position].first << "," << grid[s->A()->position].second << ");("  << grid[s->B()->position].first << "," << grid[s->B()->position].second << ")]" << std::endl;
        }

        template<typename spring_system_type>
        void set_resting_lengths(spring_system_type &spring_system, value_type time)
        {
            typedef typename spring_system_type::spring_iterator iterator;
            grid_type &grid = this->grid();
//             std::cout.precision(16);
//             std::cout << "resting lengths before = ";
//             for (iterator s = spring_system.springs_begin(), end = spring_system.springs_end(); s != end; ++s)
//             {
//                 std::cout << s->resting_length() << ",";
//             }
//             std::cout << std::endl;
            
            for (iterator s = spring_system.springs_begin(), end = spring_system.springs_end(); s != end; ++s)
            {
                size_t Ai = grid[s->A()->position].first;
                size_t Aj = grid[s->A()->position].second;
                size_t Bi = grid[s->B()->position].first;
                size_t Bj = grid[s->B()->position].second;
                s->resting_length() = m_swimmers->get_distance(Ai,Aj,Bi,Bj,time);
            }
//             std::cout << "resting lengths after = ";
//             for (iterator s = spring_system.springs_begin(), end = spring_system.springs_end(); s != end; ++s)
//             {
//                 std::cout << s->resting_length() << ",";
//             }
//             std::cout << std::endl;

        }

        template<typename vector_type>
        void get_squared_patch(vector_type &patch, size_t M, size_t N, value_type *dx)
        {
            for (size_t i = 0; i < M; ++i)
                for (size_t j = 0; j < N; ++j)
                {
                    patch.push_back(i*dx[0]);
                    patch.push_back(j*dx[1]);
                }
        }

        template<typename vector_type>
        void get_sphere_patch(vector_type &patch, size_t M, size_t N, value_type *dx)
        {
            for (size_t i = 0; i < M; ++i)
                for (size_t j = 0; j < N; ++j)
                    for (size_t l = 0; l < N; ++l)
                    {
//                     patch.push_back(i*dx[0]);
//                     patch.push_back(j*dx[1]);
//                     patch.push_back(l*dx[2]);
                    }
        }
};


template<typename _value_type>
struct surface_traits<Swimmers<_value_type> >
{
    typedef _value_type value_type;
    typedef SineGeometry<value_type> geometry_type;
};

#endif
