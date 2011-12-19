#ifndef BASE_GEOMETRY_HPP
#define BASE_GEOMETRY_HPP

#include <map>
#include<cmath>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

template<typename T> struct geometry_traits;

template<typename Derived>
class BaseGeometry
{
    protected:
        typedef typename geometry_traits<Derived>::value_type value_type;

    protected:

        vtkSmartPointer<vtkCellArray> m_cells;


    public:
        BaseGeometry() : m_cells(vtkSmartPointer<vtkCellArray>::New()) {}
        ~BaseGeometry() {}

        inline Derived *derived()
        {
            return static_cast<Derived *>(this);
        }

        template<typename array_type>
        inline void add_plane_connections(size_t i, size_t j, size_t M, size_t N, array_type &col_idx, size_t offset = 0)
        {
            int connections[8][2] = {{i - 1, j - 1},
                {i, j - 1},
                {i + 1, j - 1},
                {i + 1, j},
                {i - 1, j},
                {i + 1, j + 1},
                {i, j + 1},
                {i - 1, j + 1}
            };
            for(int k = 0; k < 8; ++k)
                if((connections[k][0] >= 0 && connections[k][1] >= 0) && (connections[k][0] <= M - 1 && connections[k][1] <= N - 1))
                    col_idx.push_back(connections[k][1]*M + connections[k][0] + offset);

        }

        template<typename array_type>
        inline void add_cylinder_connections(size_t i, size_t j, size_t M, size_t N, array_type &col_idx, size_t offset = 0)
        {
            int connections[6][2] = {{i + 1, j - 1},
                {i + 1, j},
                {i + 1, j + 1},
                {i - 1, j - 1},
                {i - 1, j},
                {i - 1, j + 1}
            };
            for(int k = 0; k < 3; ++k)
            {
                if((connections[k][1] >= 0 && connections[k][1] <= N - 1) && i == M - 1)
                    col_idx.push_back(connections[k][1]*M + (connections[k][0] + M) % M + offset);
            }
            for(int k = 3; k < 6; ++k)
            {
                if((connections[k][1] >= 0 && connections[k][1] <= N - 1) && i == 0)
                    col_idx.push_back(connections[k][1]*M + (connections[k][0] + M) % M + offset);
            }
        }


        template<typename array_type>
        inline void add_closed_connections(size_t i, size_t j, size_t M, size_t N, array_type &col_idx, size_t offset = 0)
        {

            int connections[8][2] = {{i - 1, j - 1},
                {i, j - 1},
                {i, j - 2},
                {i + 1, j - 1},
                {i + 1, j + 1},
                {i, j + 1},
                {i, j + 2},
                {i - 1, j + 1}
            };
            for(int k = 0; k < 8; ++k)
                if((connections[0][k] >= 0) && (connections[0][k] <= M - 1) && (j == 0 || j == N - 1))
                    col_idx.push_back((connections[1][k] + N) % N * M + connections[0][k] + offset);

        }

        void set_plane_cells(size_t i, size_t j, size_t M, size_t N, size_t offset = 0)
        {
            if(i < M - 1 && j < N - 1)
            {
                vtkIdType cell[4] = {j *M + i + offset, j *M + i + 1 + offset, (j + 1) *M + i + 1 + offset, (j + 1) *M + i + offset};
                m_cells->InsertNextCell(4, cell);
            }
        }

        void set_corner_cells(size_t i, size_t j, size_t M, size_t N, size_t offset = 0)
        {
            if(i == M - 1 && j < N - 1)
            {
                vtkIdType cell[4] = {j *M + i + offset, j *M + offset, (j + 1) *M + offset, (j + 1) *M + i + offset};
                m_cells->InsertNextCell(4, cell);
            }
        }

        void set_top_cells(size_t i, size_t j, size_t M, size_t N, size_t offset = 0)
        {

            if(i < M - 1 && j == N - 1)
            {
                vtkIdType cell[4] = {j *M + i + offset, j *M + i + 1 + offset, i + 1 + offset, i + offset};
                m_cells->InsertNextCell(4, cell);
            }
        }

        inline value_type get_distance(size_t Ai, size_t Aj, size_t Bi, size_t Bj, value_type scale)
        {
            value_type points[2][3] = {0};
            size_t *dims = derived()->get_dimensions();
            value_type dtheta = 2 * M_PI / dims[0];
            value_type dalpha = 2 * M_PI / dims[1];
            derived()->surface_point(Ai, Aj, scale, points[0], dtheta, dalpha);
            derived()->surface_point(Bi, Bj, scale, points[1], dtheta, dalpha);
            value_type dx[3] = {points[0][0] - points[1][0],
                                points[0][1] - points[1][1],
                                points[0][2] - points[1][2]
                               };
            return std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
        }

        template<typename spring_type>
        inline void setRestingLength(spring_type *spring, value_type time)
        {
            value_type points[2][3] = {0};
            size_t *dims = derived()->get_dimensions();
            value_type dtheta = 2 * M_PI / dims[0];
            value_type dalpha = 2 * M_PI / dims[1];
            derived()->surface_point(spring->A()->i, spring->A()->j, time, points[0], dtheta, dalpha);
            derived()->surface_point(spring->B()->i, spring->B()->j, time, points[1], dtheta, dalpha);
            value_type dx[3] = {points[0][0] - points[1][0],
            points[0][1] - points[1][1],
            points[0][2] - points[1][2]
            };
            spring->resting_length() = std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
        }

        template<typename particle_type, typename grid_type>
        void init(particle_type *particles, grid_type &grid)
        {
            size_t *dims = derived()->get_dimensions();
            value_type dtheta = 2 * M_PI / dims[0];
            value_type dalpha = 2 * M_PI / dims[1];

            for(size_t i = 0, idx = 0; i < dims[1]; ++i)
                for(size_t j = 0; j < dims[0]; ++j, ++idx)
                {
                    derived()->surface_point(i, j, 1.0, particles[idx], dtheta, dalpha);
                    grid[&particles[idx]] = std::make_pair(i, j);
                }
        }

        template<typename particle_type>
        void init(particle_type *particles, size_t num_rings)
        {
            size_t *dims = derived()->get_dimensions();
            value_type dtheta = 2 * M_PI / dims[0];
            value_type dalpha = 2 * M_PI / dims[1];

            for(size_t l = 0, idx = 0; l < num_rings; ++l)
                for(size_t i = 0; i < dims[1]; ++i)
                    for(size_t j = 0; j < dims[0]; ++j, ++idx)
                    {
                        value_type scale = .1 * l + .1;
                        derived()->surface_point(i, j, scale, particles[idx], dtheta, dalpha);
                    }
        }

        template<typename particle_type>
        void set_positions(particle_type *particles, value_type t)
        {
            size_t *dims = derived()->get_dimensions();
            value_type dtheta = 2 * M_PI / dims[0];
            value_type dalpha = 2 * M_PI / dims[1];
            for(size_t i = 0, idx = 0; i < dims[1]; ++i)
                for(size_t j = 0; j < dims[0]; ++j, ++idx)
                    derived()->surface_point(i, j, t, particles[idx], dtheta, dalpha);
        }

        void setCells()
        {
            derived()->setCells();
        }

        vtkSmartPointer<vtkCellArray> &getCells()
        {
            return m_cells;
        }

        inline void get_forcing_range(size_t &hi, size_t &lo) { derived()->get_forcing_range(hi, lo); }
        inline size_t size() { return derived()->size(); }
        inline void set_positions(value_type *positions) { derived()->set_positions(positions); }
        inline size_t const *get_dimensions() const { return derived()->get_dimensions(); }
        inline size_t *get_dimensions() { return derived()->get_dimensions(); }
};


#endif
