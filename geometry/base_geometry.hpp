#ifndef BASE_GEOMETRY_HPP
#define BASE_GEOMETRY_HPP

#include <list>
#include <map>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

template<typename T> struct geometry_traits;


class abstract_base_geometry
{
    
};

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
            return static_cast<Derived*>(this);
        }

        template<typename array_type>
        inline void add_cylinder_connections(size_t i, size_t j, array_type &col_idx)
        {
            size_t *dims = derived()->get_dimensions();
            if (i == 0 && j == 0)
            {
                col_idx.push_back(dims[0]-1);
                col_idx.push_back((j+1)*dims[0]+dims[0]-1);
            }
            if (i == dims[0]-1 && j == 0)
            {
                col_idx.push_back(j*dims[0]);
                col_idx.push_back((j+1)*dims[0]);
            }
            if (i == dims[0]-1 && j == dims[1]-1)
            {
                col_idx.push_back((j-1)*dims[0]);
                col_idx.push_back(j*dims[0]);
            }

            if (i == 0 && j == dims[1]-1)
            {
                col_idx.push_back((j-1)*dims[0]+dims[0]-1);
                col_idx.push_back(j*dims[0]+dims[0]-1);
            }
            if (i == 0 && (j > 0 && j < dims[1]-1))
            {
                col_idx.push_back((j-1)*dims[0]+dims[0]-1);
                col_idx.push_back(j*dims[0]+dims[0]-1);
                col_idx.push_back((j+1)*dims[0]+dims[0]-1);
            }
            if (i == dims[0]-1 && (j > 0 && j < dims[1]-1))
            {
                col_idx.push_back((j-1)*dims[0]);
                col_idx.push_back(j*dims[0]);
                col_idx.push_back(j*dims[0]+i+1);
            }
        }

        template<typename array_type>
        inline void add_plane_connections(size_t i, size_t j, array_type &col_idx)
        {
            size_t *dims = derived()->get_dimensions();
            if (i == 0 && j == 0)
            {
                col_idx.push_back(j*dims[0]+i+1);
                col_idx.push_back((j+1)*dims[0]+i+1);
                col_idx.push_back((j+1)*dims[0]+i);
            }
            if (i == dims[0]-1 && j == 0)
            {
                col_idx.push_back((j+1)*dims[0]+i);
                col_idx.push_back(j*dims[0]+i-1);
                col_idx.push_back((j+1)*dims[0]+i-1);
            }
            if (i == dims[0]-1 && j == dims[1]-1)
            {
                col_idx.push_back((j-1)*dims[0]+dims[0]-1);
                col_idx.push_back((j-1)*dims[0]+dims[0]-2);
                col_idx.push_back(j*dims[0]+dims[0]-2);
            }

            if (i == 0 && j == dims[1]-1)
            {
                col_idx.push_back((j-1)*dims[0]+i+1);
                col_idx.push_back(j*dims[0]+i+1);
                col_idx.push_back((j-1)*dims[0]+i);
            }
            if (i == 0 && (j > 0 && j < dims[1]-1))
            {
                col_idx.push_back((j-1)*dims[0]+i+1);
                col_idx.push_back(j*dims[0]+i+1);
                col_idx.push_back((j+1)*dims[0]+i+1);
                col_idx.push_back((j-1)*dims[0]+i);
                col_idx.push_back((j+1)*dims[0]+i);
            }
            if (i == dims[0]-1 && (j > 0 && j < dims[1]-1))
            {
                col_idx.push_back((j-1)*dims[0]+i);
                col_idx.push_back((j+1)*dims[0]+i);
                col_idx.push_back((j-1)*dims[0]+i-1);
                col_idx.push_back(j*dims[0]+i-1);
                col_idx.push_back((j+1)*dims[0]+i-1);
            }
            if (j == 0 && (i > 0 && i < dims[0]-1))
            {
                col_idx.push_back(j*dims[0]+i+1);
                col_idx.push_back((j+1)*dims[0]+i+1);
                col_idx.push_back((j+1)*dims[0]+i);
                col_idx.push_back(j*dims[0]+i-1);
                col_idx.push_back((j+1)*dims[0]+i-1);
            }
            if (j == dims[1]-1 && (i > 0 && i < dims[0]-1))
            {
                col_idx.push_back((j-1)*dims[0]+i+1);
                col_idx.push_back(j*dims[0]+i+1);
                col_idx.push_back((j-1)*dims[0]+i);
                col_idx.push_back((j-1)*dims[0]+i-1);
                col_idx.push_back(j*dims[0]+i-1);
            }
            if (i > 0 && j > 0 && (i < dims[0]-1 && j < dims[1]-1))
            {
                col_idx.push_back((j-1)*dims[0]+i+1);
                col_idx.push_back(j*dims[0]+i+1);
                col_idx.push_back((j+1)*dims[0]+i+1);
                col_idx.push_back((j-1)*dims[0]+i);
                col_idx.push_back((j+1)*dims[0]+i);
                col_idx.push_back((j-1)*dims[0]+i-1);
                col_idx.push_back(j*dims[0]+i-1);
                col_idx.push_back((j+1)*dims[0]+i-1);
            }
        }

        template<typename array_type>
        inline void add_closed_connections(size_t i, size_t j, array_type &col_idx)
        {
            size_t *dims = derived()->get_dimensions();
            size_t M = dims[0];
            size_t N = dims[1];

            if (i == 0 && j == 0)
            {
                col_idx.push_back((N-1)*M+M-1);
                col_idx.push_back((N-1)*M+i);
                col_idx.push_back((N-1)*M+i+1);
            }
            if (i == M-1 && j == 0)
            {
                col_idx.push_back((N-1)*M+M-2);
                col_idx.push_back((N-1)*M+M-1);
                col_idx.push_back((N-1)*M);
            }

            if (i == M-1 && j == N-1)
            {
                col_idx.push_back(M-2);
                col_idx.push_back(M-1);
                col_idx.push_back(0);
            }
//
            if (i == 0 && j == N-1)
            {
                col_idx.push_back(M-1);
                col_idx.push_back(0);
                col_idx.push_back(1);
            }

            if (j == 0 && (i > 0 && i < M-1))
            {
                col_idx.push_back((N-1)*M+i-1);
                col_idx.push_back((N-1)*M+i);
                col_idx.push_back((N-1)*M+i+1);
            }

            if (j == N-1 && (i > 0 && i < M-1))
            {
                col_idx.push_back(i-1);
                col_idx.push_back(i);
                col_idx.push_back(i+1);
            }

        }

        void set_plane_cells(size_t i, size_t j)
        {
            size_t *dims = derived()->get_dimensions();
            size_t M = dims[0];
            size_t N = dims[1];
            if (i < M-1 && j < N-1)
            {
                vtkIdType cell[4] = {j*M+i,j*M+i+1,(j+1)*M+i+1,(j+1)*M+i};
                m_cells->InsertNextCell(4, cell);
            }
        }

        void set_corner_cells(size_t i, size_t j)
        {
            size_t *dims = derived()->get_dimensions();
            size_t M = dims[0];
            size_t N = dims[1];
            if (i == M-1 && j < N-1)
            {
                vtkIdType cell[4] = {j*M+i,j*M,(j+1)*M,(j+1)*M+i};
                m_cells->InsertNextCell(4, cell);
            }
        }
        
        void set_top_cells(size_t i, size_t j)
        {
            size_t *dims = derived()->get_dimensions();
            size_t M = dims[0];
            size_t N = dims[1];
            if (i < M-1 && j == N-1)
            {
                vtkIdType cell[4] = {j*M+i,j*M+i+1,i+1,i};
                m_cells->InsertNextCell(4, cell);
            }
        }
        
        inline value_type get_distance(size_t Ai, size_t Aj, size_t Bi, size_t Bj, value_type t)
        {
            value_type point1[3] = {0}, point2[3] = {0};
            size_t *dims = derived()->get_dimensions();
            value_type dtheta = 2*M_PI/dims[0];
            value_type dalpha = 2*M_PI/dims[1];
            derived()->surface_point(Ai,Aj,t,point1,dtheta,dalpha);
            derived()->surface_point(Bi,Bj,t,point2,dtheta,dalpha);
            value_type dx[3] = {point1[0]-point2[0],point1[1]-point2[1],point1[2]-point2[2]};
            return std::sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
        }

        template<typename array_type>
        inline void get_connections(array_type &col_ptr, array_type &col_idx)
        {
            derived()->get_connections(col_ptr,col_idx);
        }

        template<typename grid_type>
        void init(value_type *positions, grid_type &grid)
        {
            size_t *dims = derived()->get_dimensions();
            value_type dtheta = 2*M_PI/dims[0];
            value_type dalpha = 2*M_PI/dims[1];
            
            for (size_t i = 0, idx = 0; i < dims[1]; ++i)
                for (size_t j = 0; j < dims[0]; ++j, idx+=3)
                {
                    derived()->surface_point(i,j,&positions[idx],1.0,dtheta,dalpha);
                    grid[&positions[idx]] = std::make_pair(i,j);
                }
        }

        void init(value_type *positions, size_t num_rings)
        {
            size_t *dims = derived()->get_dimensions();            
            value_type dtheta = 2*M_PI/dims[0];
            value_type dalpha = 2*M_PI/dims[1];
            
            for (size_t l = 0, idx = 0; l < num_rings; ++l)
                for (size_t i = 0; i < dims[1]; ++i)
                    for (size_t j = 0; j < dims[0]; ++j, idx+=3)
                    {
                        value_type depth = .1*l+.1;
                        derived()->surface_point(i,j,&positions[idx],depth,dtheta,dalpha);
                    }
        }

        void set_positions(value_type *positions, value_type t)
        {
            size_t *dims = derived()->get_dimensions();
            value_type dtheta = 2*M_PI/dims[0];
            value_type dalpha = 2*M_PI/dims[1];
            for (size_t i = 0, idx = 0; i < dims[1]; ++i)
                for (size_t j = 0; j < dims[0]; ++j, idx+=3)
                    derived()->surface_point(i,j,t,&positions[idx],dtheta,dalpha);
        }

        inline void setCells()
        {
            derived()->setCells();
        }

        inline vtkSmartPointer<vtkCellArray> &getCells()
        {
            return m_cells;
        }
        
        inline void get_forcing_range(size_t &hi, size_t &lo){ derived()->get_forcing_range(hi,lo); }
        inline size_t size() { return derived()->size(); }
        inline void set_positions(value_type *positions) { derived()->set_positions(positions); }
        inline size_t const *get_dimensions() const { return derived()->get_dimensions(); }
        inline size_t *get_dimensions() { return derived()->get_dimensions(); }
        inline value_type *get_x0() { return derived()->get_x0(); }
};


#endif
