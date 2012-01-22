#ifndef VTK_PARTICLE_SYSTEM_STORAGE_HPP
#define VTK_PARTICLE_SYSTEM_STORAGE_HPP
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
#include<vtkSmartPointer.h>
#include<vtkPolyData.h>
#include<vtkPointData.h>
#include<vtkPoints.h>
#include<vtkCellArray.h>

#include "particle_system_storage.hpp"

/// @name ParticleSystemStorage - Provides the storage for the particle system.
///
/// @section Description This particle system provides allocation and handling of 
/// 			the simulation data.  The purpose of this class is to provide 
/// 			storage model that is easy to swap when needing.
/// 			


template<typename vtk_data_array_type>
struct vtk_arrays
{
    vtkSmartPointer<vtk_data_array_type> positions;
    vtkSmartPointer<vtk_data_array_type> velocities;
    vtkSmartPointer<vtk_data_array_type> forces;
};


/** \internal
 * 
 * \class particle_system_storage
 *
 * \brief Stores the data of the particle system
 *
 * This class stores the data
 *
 */
template<typename T, typename particle_type, int immerse_structure_type, typename vtk_data_array_type = vtkDoubleArray>
class vtkParticleSystemStorage;

template<typename T, typename particle_type, typename vtk_data_array_type >
class vtkParticleSystemStorage<T,particle_type,PSYS::SURFACE,vtk_data_array_type>
{
    private:
        particle_system_arrays<T*,particle_type*> m_data;
        vtk_arrays<vtk_data_array_type>         m_vtk_data;
        vtkSmartPointer<vtkPolyData>            m_poly_data; 

    public:
        inline explicit vtkParticleSystemStorage(size_t data_size) : m_poly_data(vtkSmartPointer<vtkPolyData>::New())
        {
            size_t num_particles = data_size/3;
            m_data.positions = new T[data_size];
            m_data.velocities = new T[data_size];
            m_data.forces = new T[data_size];
            std::fill(m_data.positions,m_data.positions+data_size,T(0));
            std::fill(m_data.velocities,m_data.velocities+data_size,T(0));
            std::fill(m_data.forces,m_data.forces+data_size,T(0));
            m_data.particles = new particle_type[num_particles];
            for (size_t i = 0, idx = 0; i < num_particles; ++i, idx+=3)
            {
                m_data.particles[i].position = &m_data.positions[idx];
                m_data.particles[i].velocity = &m_data.velocities[idx];
                m_data.particles[i].force = &m_data.forces[idx];
            }

            m_vtk_data.positions = vtkSmartPointer<vtk_data_array_type>::New();
            m_vtk_data.velocities = vtkSmartPointer<vtk_data_array_type>::New();
            m_vtk_data.forces = vtkSmartPointer<vtk_data_array_type>::New();
            
            m_vtk_data.positions->SetArray(m_data.positions,data_size,1);
            m_vtk_data.velocities->SetArray(m_data.velocities,data_size,1);
            m_vtk_data.forces->SetArray(m_data.forces,data_size,1);
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

        ~vtkParticleSystemStorage()
        {
            delete [] m_data.positions;
            delete [] m_data.velocities;
            delete [] m_data.forces;
            delete [] m_data.particles;
        }
        inline void swap(vtkParticleSystemStorage& other) { std::swap(m_data,other.m_data); }

        inline const T *position(size_t i) const { return m_data.particle[i].position; }
        inline T *position(size_t i)             { return m_data.particle[i].position; }

        inline const T *velocity(size_t i) const { return m_data.particle[i].velocity; }
        inline T *velocity(size_t i)             { return m_data.particle[i].velocity; }

        inline const T *force(size_t i) const { return m_data.particle[i].force; }
        inline T *force(size_t i)             { return m_data.particle[i].force; }

        inline const particle_type *particles(size_t i) const { return m_data.particles[i]; }
        inline particle_type *particles(size_t i)             { return m_data.particles[i]; }

        inline const T *positions() const { return m_data.positions; }
        inline T *positions()             { return m_data.positions; }

        inline const T *velocities() const { return m_data.velocities; }
        inline T *velocities()             { return m_data.velocities; }

        inline const T *forces() const { return m_data.forces; }
        inline T *forces()             { return m_data.forces; }

        inline const particle_type *particles() const { return m_data.particles; }
        inline particle_type *particles()             { return m_data.particles; }

        
        inline vtkSmartPointer<vtkPolyData> &grid() { return m_poly_data; }
        
};

template<typename T, typename particle_type, typename vtk_data_array_type>
class vtkParticleSystemStorage<T,particle_type,PSYS::VOLUME,vtk_data_array_type>
{
    private:
        particle_system_arrays<T*,particle_type*> m_data;
        vtk_arrays<vtk_data_array_type>         m_vtk_data;
        vtkSmartPointer<vtkPolyData>            m_poly_data; 

    public:
        inline explicit vtkParticleSystemStorage(size_t data_size) : m_poly_data(vtkSmartPointer<vtkPolyData>::New())
        {
            size_t num_particles = data_size/3;
            m_data.positions = new T[data_size];
            m_data.velocities = new T[data_size];
            m_data.forces = 0;
            m_data.particles = new particle_type[num_particles];
            for (size_t i = 0, idx = 0; i < num_particles; ++i, idx+=3)
            {
                m_data.particles[i].position = &m_data.positions[idx];
                m_data.particles[i].velocity = &m_data.velocities[idx];
            }

            m_vtk_data.positions = vtkSmartPointer<vtk_data_array_type>::New();
            m_vtk_data.velocities = vtkSmartPointer<vtk_data_array_type>::New();
            m_vtk_data.forces = vtkSmartPointer<vtk_data_array_type>::New();
            
            m_vtk_data.positions->SetArray(m_data.positions,data_size,1);
            m_vtk_data.velocities->SetArray(m_data.velocities,data_size,1);
            m_vtk_data.positions->SetNumberOfComponents(3);
            m_vtk_data.positions->SetName("positions");
            m_vtk_data.velocities->SetNumberOfComponents(3);
            m_vtk_data.velocities->SetName("velocity");

            vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
            points->SetData(m_vtk_data.positions);
            m_poly_data->SetPoints(points);
            
            if (m_vtk_data.velocities->GetSize() > 0)
                m_poly_data->GetPointData()->AddArray(m_vtk_data.velocities);

            vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();            
            for (vtkIdType i = 0; i < num_particles; ++i)
            {
                cells->InsertNextCell(VTK_VERTEX,&i);
            }
            m_poly_data->SetVerts(cells);
        }

        ~vtkParticleSystemStorage()
        {
            delete [] m_data.positions;
            delete [] m_data.velocities;
            delete [] m_data.particles;
        }
        inline void swap(vtkParticleSystemStorage& other) { std::swap(m_data,other.m_data); }

        inline const T *position(size_t i) const { return m_data.particle[i].position; }
        inline T *position(size_t i)             { return m_data.particle[i].position; }

        inline const T *velocity(size_t i) const { return m_data.particle[i].velocity; }
        inline T *velocity(size_t i)             { return m_data.particle[i].velocity; }

        inline const particle_type *particles(size_t i) const { return m_data.particles[i]; }
        inline particle_type *particles(size_t i)             { return m_data.particles[i]; }

        inline const T *force(size_t i) const { return m_data.particle[i].force; }
        inline T *force(size_t i)             { return m_data.particle[i].force; }

        inline const T *positions() const { return m_data.positions; }
        inline T *positions()             { return m_data.positions; }

        inline const T *forces() const { return m_data.forces; }
        inline T *forces()             { return m_data.forces; }

        inline const T *velocities() const { return m_data.velocities; }
        inline T *velocities()             { return m_data.velocities; }

        inline const particle_type *particles() const { return m_data.particles; }
        inline particle_type *particles()             { return m_data.particles; }

        inline vtkSmartPointer<vtkPolyData> &grid() { return m_poly_data; }
};


#endif
