#ifndef VTK_STORAGE_WRAPPER_HPP
#define VTK_STORAGE_WRAPPER_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    vtkParticleSystemStorage
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================

#include<vtkDataArray.h>
#include<vtkDoubleArray.h>
#include<vtkFloatArray.h>
#include<vtkSmartPointer.h>
#include<vtkPolyData.h>
#include<vtkPointData.h>
#include<vtkPoints.h>
#include<vtkCellArray.h>

template<typename vtk_data_array_type>
struct vtkArrays
{
    vtkSmartPointer<vtk_data_array_type> positions;
    vtkSmartPointer<vtk_data_array_type> velocities;
    vtkSmartPointer<vtk_data_array_type> forces;
    vtkSmartPointer<vtkCellArray> cells;
};


template<typename particle_system_storage, typename array_type = vtkDoubleArray >
class vtkStorageWrapper
{
    private:
        vtkArrays<array_type>         m_vtk_data;
        vtkSmartPointer<vtkPolyData>  m_poly_data;
        vtkIdType                     m_id_buffer[4][2];

    public:
        vtkStorageWrapper(particle_system_storage &data) : m_poly_data(vtkSmartPointer<vtkPolyData>::New())
        {
            m_vtk_data.positions = vtkSmartPointer<array_type>::New();
            m_vtk_data.velocities = vtkSmartPointer<array_type>::New();
            m_vtk_data.forces = vtkSmartPointer<array_type>::New();
            m_vtk_data.cells = vtkSmartPointer<vtkCellArray>::New();
            
            size_t size = data.data_size();
            m_vtk_data.positions->SetArray(data.positions(), size, 1);
            m_vtk_data.velocities->SetArray(data.velocities(), size, 1);
            m_vtk_data.forces->SetArray(data.forces(), size, 1);
            m_vtk_data.positions->SetNumberOfComponents(3);
            m_vtk_data.positions->SetName("positions");
            m_vtk_data.velocities->SetNumberOfComponents(3);
            m_vtk_data.velocities->SetName("velocity");
            m_vtk_data.forces->SetNumberOfComponents(3);
            m_vtk_data.forces->SetName("force");

            vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
            points->SetData(m_vtk_data.positions);
            m_poly_data->SetPoints(points);
            if (m_vtk_data.velocities->GetSize() > 0)
                m_poly_data->GetPointData()->AddArray(m_vtk_data.velocities);
            if (m_vtk_data.forces->GetSize() > 0)
                m_poly_data->GetPointData()->AddArray(m_vtk_data.forces);
            m_poly_data->SetPolys(m_vtk_data.cells);
        }

        inline const vtkSmartPointer<array_type> &positions() const    { return m_vtk_data.positions; }
        inline vtkSmartPointer<array_type> &positions()                { return m_vtk_data.positions; }
        inline const vtkSmartPointer<array_type> &velocities() const   { return m_vtk_data.velocities; }
        inline vtkSmartPointer<array_type> &velocities()               { return m_vtk_data.velocities; }
        inline const vtkSmartPointer<array_type> &forces() const       { return m_vtk_data.forces; }
        inline vtkSmartPointer<array_type> &forces()                   { return m_vtk_data.forces; }
        inline const vtkSmartPointer<vtkCellArray> &cells() const      { return m_vtk_data.cells; }
        inline vtkSmartPointer<vtkCellArray> &cells()                  { return m_vtk_data.cells; }
        inline const vtkSmartPointer<vtkPolyData> &grid() const        { return m_poly_data; }
        inline vtkSmartPointer<vtkPolyData> &grid()                    { return m_poly_data; }

        inline void setInnerCells(int i, int j, int M, int N, size_t offset)
        {
            if(i < M - 1 && j < N - 1)
            {
                setIdBuffer(i, j);
                vtkIdType cell[4] =
                {
                    m_id_buffer[0][1]*M + m_id_buffer[0][0] + offset,
                    m_id_buffer[1][1]*M + m_id_buffer[1][0] + offset,
                    m_id_buffer[2][1]*M + m_id_buffer[2][0] + offset,
                    m_id_buffer[3][1]*M + m_id_buffer[3][0] + offset
                };
                m_vtk_data.cells->InsertNextCell(4, cell);
            }
        }

        inline void setEdgeCells(int i, int j, int M, int N, size_t offset)
        {
            if(i == M - 1 && j < N - 1)
            {
                setIdBuffer(i, j);
                vtkIdType cell[4] =
                {
                    m_id_buffer[0][1]*M + m_id_buffer[0][0] + offset,
                    m_id_buffer[1][1]*M + m_id_buffer[1][0]%M + offset,
                    m_id_buffer[2][1]*M + m_id_buffer[2][0]%M + offset,
                    m_id_buffer[3][1]*M + m_id_buffer[3][0] + offset
                };
                m_vtk_data.cells->InsertNextCell(4, cell);
            }
        }

        inline void setTopCells(int i, int j, int M, int N, size_t offset)
        {
            if(i <= M - 1 && j == N - 1)
            {
                setIdBuffer(i, j);
                vtkIdType cell[4] =
                {
                    m_id_buffer[0][1]*M + m_id_buffer[0][0] + offset,
                    m_id_buffer[1][1]*M + m_id_buffer[1][0]%M + offset,
                    m_id_buffer[2][1]%N*M + m_id_buffer[2][0]%M + offset,
                    m_id_buffer[3][1]%N*M + m_id_buffer[3][0] + offset
                };
                m_vtk_data.cells->InsertNextCell(4, cell);
            }
        }

        void setEdgeCells(int M, int N, size_t offset = 0)
        {
            for(int j = 0; j < N; ++j)
                for(int i = 0; i < M; ++i)
                    setEdgeCells(i,j,M,N,offset);
        }

        void setTopCells(int M, int N, size_t offset = 0)
        {
            for(int j = 0; j < N; ++j)
                for(int i = 0; i < M; ++i)
                    setTopCells(i,j,M,N,offset);
        }

        void setInnerCells(int M, int N, size_t offset = 0)
        {
            for(int j = 0; j < N; ++j)
                for(int i = 0; i < M; ++i)
                    setInnerCells(i,j,M,N,offset);
        }
    private:
        inline void setIdBuffer(int i, int j)
        {
            m_id_buffer[0][0] = i;   m_id_buffer[0][1]   = j;
            m_id_buffer[1][0] = i + 1; m_id_buffer[1][1] = j;
            m_id_buffer[2][0] = i + 1; m_id_buffer[2][1] = j + 1;
            m_id_buffer[3][0] = i;   m_id_buffer[3][1]   = j + 1;
        }

};

#endif
