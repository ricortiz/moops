#ifndef BASE_GEOMETRY_HPP
#define BASE_GEOMETRY_HPP

#include <map>
#include<cmath>
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
        inline void add_plane_connections(size_t i, size_t j, size_t M, size_t N, array_type &col_idx, size_t offset = 0)
        {

            int connections[10][2] = {{i-1,j-1},
                                      {i,j-1},
                                      {i+1,j-1},
                                      {i+2,j},
                                      {i+1,j},
                                      {i-1,j},
                                      {i-2,j},
                                      {i+1,j+1},
                                      {i,j+1},
                                      {i-1,j+1}};
            for(int k = 0; k < 10; ++k)
                if((connections[k][0] >=0 && connections[k][1] >=0) && (connections[k][0] <= M-1 && connections[k][1] <= N-1))
                    col_idx.push_back(connections[k][1]*M+connections[k][0]+offset);
//             if (i == 0 && j == 0)
//             {
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back((j+1)*M+i);
//                 col_idx.push_back((j+1)*M+i+1);
//             }
//             if (i == M-1 && j == 0)
//             {
//                 col_idx.push_back((j+1)*M+i);
//                 col_idx.push_back(j*M+i-1);
//                 col_idx.push_back((j+1)*M+i-1);
//             }
//             if (i == M-1 && j == N-1)
//             {
//                 col_idx.push_back((j-1)*M+M-1);
//                 col_idx.push_back((j-1)*M+M-2);
//                 col_idx.push_back(j*M+M-2);
//             }
//             
//             if (i == 0 && j == N-1)
//             {
//                 col_idx.push_back((j-1)*M+i+1);
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back((j-1)*M+i);
//             }
//             if (i == 0 && (j > 0 && j < N-1))
//             {
//                 col_idx.push_back((j-1)*M+i+1);
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back((j+1)*M+i+1);
//                 col_idx.push_back((j-1)*M+i);
//                 col_idx.push_back((j+1)*M+i);
//             }
//             if (i == M-1 && (j > 0 && j < N-1))
//             {
//                 col_idx.push_back((j-1)*M+i);
//                 col_idx.push_back((j+1)*M+i);
//                 col_idx.push_back((j-1)*M+i-1);
//                 col_idx.push_back(j*M+i-1);
//                 col_idx.push_back((j+1)*M+i-1);
//             }
//             if (j == 0 && (i > 0 && i < M-1))
//             {
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back((j+1)*M+i+1);
//                 col_idx.push_back((j+1)*M+i);
//                 col_idx.push_back(j*M+i-1);
//                 col_idx.push_back((j+1)*M+i-1);
//             }
//             if (j == N-1 && (i > 0 && i < M-1))
//             {
//                 col_idx.push_back((j-1)*M+i+1);
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back((j-1)*M+i);
//                 col_idx.push_back((j-1)*M+i-1);
//                 col_idx.push_back(j*M+i-1);
//             }
//             if (i > 0 && j > 0 && (i < M-1 && j < N-1))
//             {
//                 col_idx.push_back((j-1)*M+i+1);
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back((j+1)*M+i+1);
//                 col_idx.push_back((j-1)*M+i);
//                 col_idx.push_back((j+1)*M+i);
//                 col_idx.push_back((j-1)*M+i-1);
//                 col_idx.push_back(j*M+i-1);
//                 col_idx.push_back((j+1)*M+i-1);
//             }
        }
        
        template<typename array_type>
        inline void add_cylinder_connections(size_t i, size_t j, size_t M, size_t N, array_type &col_idx, size_t offset = 0)
        {
            int connections[4][2] = {{i+2,j},
                                     {i+1,j},
                                     {i-1,j},
                                     {i-2,j}};
            for(int k = 0; k < 4; ++k)
                if((connections[k][0] < 0 && connections[k][1] >= 0) || (connections[k][0] > M-1 && connections[k][1] <= N-1))
                    col_idx.push_back(connections[k][1]*M+connections[k][0]%M+offset);
//                     col_idx.push_back(j*M+(i+2)%M);
//             if (i == 0 && j == 0)
//             {
//                 col_idx.push_back(M-1);
//                 col_idx.push_back((j+1)*M+M-1);
//             }
//             if (i == M-1 && j == 0)
//             {
//                 col_idx.push_back(j*M);
//                 col_idx.push_back((j+1)*M);
//             }
//             if (i == M-1 && j == N-1)
//             {
//                 col_idx.push_back((j-1)*M);
//                 col_idx.push_back(j*M);
//             }
//             
//             if (i == 0 && j == N-1)
//             {
//                 col_idx.push_back((j-1)*M+M-1);
//                 col_idx.push_back(j*M+M-1);
//             }
//             if (i == 0 && (j > 0 && j < N-1))
//             {
//                 col_idx.push_back((j-1)*M+M-1);
//                 col_idx.push_back(j*M+M-1);
//                 col_idx.push_back((j+1)*M+M-1);
//             }
//             if (i == M-1 && (j > 0 && j < N-1))
//             {
//                 col_idx.push_back((j-1)*M);
//                 col_idx.push_back(j*M);
//                 col_idx.push_back(j*M+i+1);
//             }
        }            
        
        template<typename array_type>
        inline void add_closed_connections(size_t i, size_t j, size_t M, size_t N, array_type &col_idx)
        {
            int connections[2][6] = {{i-1,i,i+1,i+1,i,i-1},{j-1,j-1,j-1,j+1,j+1,j+1}};
            for(int k = 0; k < 6; ++k)
                if((connections[0][k] >= 0 && connections[1][k] < 0) || (connections[0][k] <= M-1 && connections[1][k] > N-1))
                    col_idx.push_back(connections[1][k]%M+connections[0][k]);
//             col_idx.push_back(((j+2)%N)*M+i);
//             col_idx.push_back(((j+3)%N)*M+i);
//             if (i == 0 && j == 0)
//             {
//                 col_idx.push_back((N-1)*M+M-1);
//                 col_idx.push_back((N-1)*M+i);
//                 col_idx.push_back((N-1)*M+i+1);
//             }
//             if (i == M-1 && j == 0)
//             {
//                 col_idx.push_back((N-1)*M+M-2);
//                 col_idx.push_back((N-1)*M+M-1);
//                 col_idx.push_back((N-1)*M);
//             }
//             
//             if (i == M-1 && j == N-1)
//             {
//                 col_idx.push_back(M-2);
//                 col_idx.push_back(M-1);
//                 col_idx.push_back(0);
//             }
//             //
//             if (i == 0 && j == N-1)
//             {
//                 col_idx.push_back(M-1);
//                 col_idx.push_back(0);
//                 col_idx.push_back(1);
//             }
//             
//             if (j == 0 && (i > 0 && i < M-1))
//             {
//                 col_idx.push_back((N-1)*M+i-1);
//                 col_idx.push_back((N-1)*M+i);
//                 col_idx.push_back((N-1)*M+i+1);
//             }
//             
//             if (j == N-1 && (i > 0 && i < M-1))
//             {
//                 col_idx.push_back(i-1);
//                 col_idx.push_back(i);
//                 col_idx.push_back(i+1);
//             }
            
        }
        
//         template<typename array_type>
//         inline void add_plane_connections(size_t i, size_t j, size_t M, size_t N, array_type &col_idx)
//         {
//             size_t offset = M/2;
//             col_idx.push_back(j*M+(i+offset)%M);
//             if (i == 0 && j == 0)
//             {
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back((j+1)*M+i+1);
//                 col_idx.push_back((j+1)*M+i);
//                 
//                 col_idx.push_back(j*M+i+2);
//                 if (N > 2)
//                     col_idx.push_back((j+2)*M+i);
//             }
//             if (i == M-1 && j == 0)
//             {
//                 col_idx.push_back((j+1)*M+i);
//                 col_idx.push_back(j*M+i-1);
//                 col_idx.push_back((j+1)*M+i-1);
//                 
//                 col_idx.push_back(j*M+i-2);
//                 if (N > 2)
//                     col_idx.push_back((j+2)*M+i);
//             }
//             if (i == M-1 && j == N-1)
//             {
//                 col_idx.push_back((j-1)*M+i);
//                 col_idx.push_back((j-1)*M+i-1);
//                 col_idx.push_back(j*M+i-1);
//                 col_idx.push_back(j*M+i-2);
//                 if (N > 2)
//                     col_idx.push_back((j-2)*M+i);
//             }
//             
//             if (i == 0 && j == N-1)
//             {
//                 col_idx.push_back((j-1)*M+i+1);
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back(j*M+i+2);
//                 col_idx.push_back((j-1)*M+i);
//                 if (N > 2)
//                     col_idx.push_back((j-2)*M+i);
//             }
//             if (i == 0 && (j > 0 && j < N-1))
//             {
//                 col_idx.push_back((j-1)*M+i+1);
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back((j+1)*M+i+1);
//                 col_idx.push_back((j-1)*M+i);
//                 col_idx.push_back((j+1)*M+i);
//                 
//                 col_idx.push_back(j*M+i+2);
//                 if( j > 1 )
//                     col_idx.push_back((j-2)*M+i);
//                 if( j < N-2 )
//                     col_idx.push_back((j+2)*M+i);
//             }
//             if (i == M-1 && (j > 0 && j < N-1))
//             {
//                 col_idx.push_back((j-1)*M+i);
//                 col_idx.push_back((j+1)*M+i);
//                 col_idx.push_back((j-1)*M+i-1);
//                 col_idx.push_back(j*M+i-1);
//                 col_idx.push_back((j+1)*M+i-1);
//                 if( j > 1 )
//                     col_idx.push_back((j-2)*M+i);
//                 if( j < N-2 )
//                     col_idx.push_back((j+2)*M+i);
//             }
//             if (j == 0 && (i > 0 && i < M-1))
//             {
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back((j+1)*M+i+1);
//                 col_idx.push_back((j+1)*M+i);
//                 col_idx.push_back(j*M+i-1);
//                 col_idx.push_back((j+1)*M+i-1);
//                 if(N > 2)
//                     col_idx.push_back((j+2)*M+i);
//                 if(i < M-2)
//                     col_idx.push_back(j*M+i+2);
//                 if(i > 1)
//                     col_idx.push_back(j*M+i-2);
//             }
//             if (j == N-1 && (i > 0 && i < M-1))
//             {
//                 col_idx.push_back((j-1)*M+i+1);
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back((j-1)*M+i);
//                 col_idx.push_back((j-1)*M+i-1);
//                 col_idx.push_back(j*M+i-1);
//                 if(N > 2)
//                     col_idx.push_back((j-2)*M+i);
//                 if(i < M-2)
//                     col_idx.push_back(j*M+i+2);
//                 if(i > 1)
//                     col_idx.push_back(j*M+i-2);
//             }
//             if (i > 0 && i < M-1 && (j > 0 && j < N-1))
//             {
//                 col_idx.push_back((j-1)*M+i+1);
//                 col_idx.push_back(j*M+i+1);
//                 col_idx.push_back((j+1)*M+i+1);
//                 col_idx.push_back((j-1)*M+i);
//                 col_idx.push_back((j+1)*M+i);
//                 col_idx.push_back((j-1)*M+i-1);
//                 col_idx.push_back(j*M+i-1);
//                 col_idx.push_back((j+1)*M+i-1);
//                 if(i < M-2)
//                     col_idx.push_back(j*M+i+2);
//                 if(i > 1)
//                     col_idx.push_back(j*M+i-2);
//                 if(j < N-2)
//                     col_idx.push_back((j+2)*M+i);
//                 if(j > 2)
//                     col_idx.push_back((j-2)*M+i);
//                     
//             }
//         }
//         
//         template<typename array_type>
//         inline void add_cylinder_connections(size_t i, size_t j,size_t M, size_t N, array_type &col_idx)
//         {
//             if (i == 0 && j == 0)
//             {
//                 col_idx.push_back(j*M+M-1);
//                 col_idx.push_back(j*M+M-2);
//                 col_idx.push_back((j+1)*M+M-1);
//             }
//             if (i == M-1 && j == 0)
//             {
//                 col_idx.push_back(j*M);
//                 col_idx.push_back(j*M+1);
//                 col_idx.push_back((j+1)*M);
//             }
//             if (i == M-1 && j == N-1)
//             {
//                 col_idx.push_back((j-1)*M);
//                 col_idx.push_back(j*M);
//                 col_idx.push_back(j*M+1);
//             }
// 
//             if (i == 0 && j == N-1)
//             {
//                 col_idx.push_back((j-1)*M+M-1);
//                 col_idx.push_back(j*M+M-1);
//                 col_idx.push_back(j*M+M-2);
//             }
//             if (i == 0 && (j > 0 && j < N-1))
//             {
//                 col_idx.push_back((j-1)*M+M-1);
//                 col_idx.push_back(j*M+M-1);
//                 col_idx.push_back(j*M+M-2);
//                 col_idx.push_back((j+1)*M+M-1);
//             }
//             if (i == M-1 && (j > 0 && j < N-1))
//             {
//                 col_idx.push_back((j-1)*M);
//                 col_idx.push_back(j*M);
//                 col_idx.push_back(j*M+1);
//                 col_idx.push_back((j+1)*M);
//             }
//         }
// 
//         template<typename array_type>
//         inline void add_closed_connections(size_t i, size_t j, size_t M, size_t N, array_type &col_idx)
//         {
//             if (i == 0 && j == 0)
//             {
//                 col_idx.push_back((N-1)*M+M-1);
//                 col_idx.push_back((N-1)*M+i);
//                 col_idx.push_back((N-1)*M+i+1);
//                 col_idx.push_back((N-2)*M+i);
//             }
//             if (i == M-1 && j == 0)
//             {
//                 col_idx.push_back((N-1)*M+M-2);
//                 col_idx.push_back((N-1)*M+M-1);
//                 col_idx.push_back((N-1)*M);
//                 col_idx.push_back((N-2)*M+M-1);
//             }
// 
//             if (i == M-1 && j == N-1)
//             {
//                 col_idx.push_back(M-2);
//                 col_idx.push_back(M-1);
//                 col_idx.push_back(0);
//                 col_idx.push_back(M+M-1);
//             }
// //
//             if (i == 0 && j == N-1)
//             {
//                 col_idx.push_back(M-1);
//                 col_idx.push_back(0);
//                 col_idx.push_back(1);
//                 col_idx.push_back(M);
//             }
// 
//             if (j == 0 && (i > 0 && i < M-1))
//             {
//                 col_idx.push_back((N-1)*M+i-1);
//                 col_idx.push_back((N-1)*M+i);
//                 col_idx.push_back((N-1)*M+i+1);
//                 col_idx.push_back((N-2)*M+i);
//             }
// 
//             if (j == N-1 && (i > 0 && i < M-1))
//             {
//                 col_idx.push_back(i-1);
//                 col_idx.push_back(i);
//                 col_idx.push_back(i+1);
//                 col_idx.push_back(M+i);
//             }
// 
//         }

        void set_plane_cells(size_t i, size_t j, size_t M, size_t N, size_t offset = 0)
        {
//             size_t *dims = derived()->get_dimensions();
//             size_t M = M;
//             size_t N = N;
            if (i < M-1 && j < N-1)
            {
                vtkIdType cell[4] = {j*M+i+offset,j*M+i+1+offset,(j+1)*M+i+1+offset,(j+1)*M+i+offset};
                m_cells->InsertNextCell(4, cell);
            }
        }

        void set_corner_cells(size_t i, size_t j, size_t M, size_t N, size_t offset = 0)
        {
//             size_t *dims = derived()->get_dimensions();
//             size_t M = M;
//             size_t N = N;
            if (i == M-1 && j < N-1)
            {
                vtkIdType cell[4] = {j*M+i+offset,j*M+offset,(j+1)*M+offset,(j+1)*M+i+offset};
                m_cells->InsertNextCell(4, cell);
            }
        }
        
        void set_top_cells(size_t i, size_t j, size_t M, size_t N, size_t offset = 0)
        {
//             size_t *dims = derived()->get_dimensions();
//             size_t M = M;
//             size_t N = N;
            if (i < M-1 && j == N-1)
            {
                vtkIdType cell[4] = {j*M+i+offset,j*M+i+1+offset,i+1+offset,i+offset};
                m_cells->InsertNextCell(4, cell);
            }
        }
        
        inline value_type get_distance(size_t Ai, size_t Aj, size_t Bi, size_t Bj, value_type scale)
        {
            value_type point1[3] = {0}, point2[3] = {0};
            size_t *dims = derived()->get_dimensions();
            value_type dtheta = 2*M_PI/dims[0];
            value_type dalpha = 2*M_PI/dims[1];
            derived()->surface_point(Ai,Aj,scale,point1,dtheta,dalpha);
            derived()->surface_point(Bi,Bj,scale,point2,dtheta,dalpha);
            value_type dx[3] = {point1[0]-point2[0],point1[1]-point2[1],point1[2]-point2[2]};
            return std::sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
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
                    derived()->surface_point(i,j,1.0,&positions[idx],dtheta,dalpha);
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
                        value_type scale = .1*l+.1;
                        derived()->surface_point(i,j,scale,&positions[idx],dtheta,dalpha);
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
