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
};


template<typename particle_system_storage, typename array_type = vtkDoubleArray >
class vtkStorageWrapper
{
    private:
        vtkArrays<array_type>         m_vtk_data;
        vtkSmartPointer<vtkPolyData>  m_poly_data;

    public:
        vtkStorageWrapper(particle_system_storage &m_data) : m_poly_data(vtkSmartPointer<vtkPolyData>::New())
        {
            m_vtk_data.positions = vtkSmartPointer<vtk_data_array_type>::New();
            m_vtk_data.velocities = vtkSmartPointer<vtk_data_array_type>::New();
            m_vtk_data.forces = vtkSmartPointer<vtk_data_array_type>::New();

            m_vtk_data.positions->SetArray(m_data.positions(), m_data.data_size(), 1);
            m_vtk_data.velocities->SetArray(m_data.velocities(), m_data.data_size(), 1);
            m_vtk_data.forces->SetArray(m_data.forces(), m_data.data_size(), 1);
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
        }

        inline const vtkSmartPointer<vtk_data_array_type> &positions() const { return m_vtk_data.positions; }
        inline vtkSmartPointer<vtk_data_array_type> &positions()             { return m_vtk_data.positions; }

        inline const vtkSmartPointer<vtk_data_array_type> &velocities() const { return m_vtk_data.velocities; }
        inline vtkSmartPointer<vtk_data_array_type> &velocities()             { return m_vtk_data.velocities; }

        inline const vtkSmartPointer<vtk_data_array_type> &forces() const { return m_vtk_data.forces; }
        inline vtkSmartPointer<vtk_data_array_type> &forces()             { return m_vtk_data.forces; }

        inline vtkSmartPointer<vtkPolyData> &grid() { return m_poly_data; }

};

#endif
