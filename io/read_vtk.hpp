#ifndef READ_VTK_H
#define READ_VTK_H
/// C++ Interfaces for reading vtk
///
/// @brief Description: Read vtk formated data
///
/// Author: Ricardo Ortiz <ricardo.ortiz@tulane.edu>, (C) 2008
/// $Id $
#include <OpenTissue/core/math/math_vector3.h>
#include <OpenTissue/core/containers/mesh/mesh.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkDataArray.h>
namespace IO
{

    template<typename particle_system>
    bool read_vtk ( int grid_form, particle_system &particles, const std::string &file_name )
    {
        typedef typename particle_system::particle_iterator particle_iterator;
        typedef typename particle_system::particle_type     particle_type;
        typedef typename particle_system::vector3_type      vector3_type;

        switch ( grid_form )
        {

            case 0:
            {
                vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
                reader->SetFileName ( file_name.c_str() );
                reader->Update();
                vtkUnstructuredGrid *vtk_grid = reader->GetOutput();
                vtk_grid->Update();

                vtkDoubleArray* velocities = static_cast<vtkDoubleArray*> ( vtk_grid->GetPointData()->GetVectors ( "velocities" ) );
                vtkDoubleArray* forces = static_cast<vtkDoubleArray*> ( vtk_grid->GetPointData()->GetVectors ( "forces" ) );

                if ( particles.particles_size() > 0 )
                {
                    particle_iterator p = particles.particle_begin();
                    particle_iterator end = particles.particle_end();

                    for ( vtkIdType i = 0;p != end, i < vtk_grid->GetNumberOfPoints();++p, ++i )
                    {
                        double * point = vtk_grid->GetPoint ( i );
                        p->position() = vector3_type ( point[0], point[1], point[2] );
                        velocities->GetTupleValue ( i, point );
                        p->velocity() = vector3_type ( point[0], point[1], point[2] );
                        forces->GetTupleValue ( i, point );
                        p->force() = vector3_type ( point[0], point[1], point[2] );
                    }
                }
                else
                {
                    for ( vtkIdType i = 0; i < vtk_grid->GetNumberOfPoints(); ++i )
                    {
                        particle_iterator p = particles.create_particle ( particle_type() );
                        double * point = vtk_grid->GetPoint ( i );
                        p->position() = vector3_type ( point[0], point[1], point[2] );
                        velocities->GetTupleValue ( i, point );
                        p->velocity() = vector3_type ( point[0], point[1], point[2] );
                        forces->GetTupleValue ( i, point );
                        p->force() = vector3_type ( point[0], point[1], point[2] );
                    }
                }

                reader->Delete();

                break;
            }

            case 1:
            {
                vtkXMLStructuredGridReader *reader = vtkXMLStructuredGridReader::New();
                reader->SetFileName ( file_name.c_str() );
                reader->Update();
                vtkStructuredGrid *vtk_grid = reader->GetOutput();
                vtk_grid->Update();

                vtkDoubleArray* velocities = static_cast<vtkDoubleArray*> ( vtk_grid->GetPointData()->GetVectors ( "velocities" ) );

                vtkDoubleArray* forces = static_cast<vtkDoubleArray*> ( vtk_grid->GetPointData()->GetVectors ( "forces" ) );

                particle_iterator p = particles.particle_begin();
                particle_iterator end = particles.particle_end();

                for ( vtkIdType i = 0;p != end, i < vtk_grid->GetNumberOfPoints();++p, ++i )
                {
                    double * point = vtk_grid->GetPoint ( i );
                    p->position() = vector3_type ( point[0], point[1], point[2] );
                    velocities->GetTupleValue ( i, point );
                    p->velocity() = vector3_type ( point[0], point[1], point[2] );
                    forces->GetTupleValue ( i, point );
                    p->force() = vector3_type ( point[0], point[1], point[2] );
                }

                reader->Delete();

                break;
            }

            default:
            {
                std::cerr << "Grid type not implemented." << std::endl;
                return 0;
            }
        }

        return 1;

    }


    namespace polymesh = OpenTissue::polymesh;

// template<>
// bool read_vtk<polymesh::PolyMesh<> >(int grid_form, polymesh::PolyMesh<> &polymesh, const std::string &file_name)
// {
//     typedef polymesh::PolyMesh<>::math_types math_types;
//     typedef math_types::vector3_type      vector3_type;
//
//     vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
//     reader->SetFileName(file_name.c_str());
//     reader->Update();
//     vtkUnstructuredGrid *vtk_grid = reader->GetOutput();
//     vtk_grid->Update();
//
//     for (vtkIdType i = 0; i < vtk_grid->GetNumberOfPoints(); ++i)
//     {
//         double * point = vtk_grid->GetPoint(i);
//         vector3_type v(point[0], point[1], point[2]);
//         polymesh.add_vertex(v);
//     }
//
//     for (vtkIdType i = 0; i < vtk_grid->GetNumberOfCells(); ++i)
//     {
//         vtkCell* cell = vtk_grid->GetCell(i);
//         int size = cell->GetNumberOfPoints();
//         vtkIdType *points = cell->GetPointIds()->GetPointer(0);
// //         std::cout << "[" << points[0] << "," << points[1] << "," << points[2] << "," << points[3] << "]" << std::endl;
//         std::list<polymesh::PolyMesh<>::vertex_handle> handles;
//
//         for (int j = 0; j < size; ++j)
//             handles.push_back(polymesh.get_vertex_handle(points[j]));
//
//         polymesh.add_face(handles.begin(), handles.end());
//     }
//
//     reader->Delete();
// }


    namespace trimesh = OpenTissue::trimesh;

// template<>
// bool read_vtk<trimesh::TriMesh<> >(int grid_form, trimesh::TriMesh<> &mesh, const std::string &file_name)
// {
//     typedef trimesh::TriMesh<>::math_types math_types;
//     typedef math_types::vector3_type      vector3_type;
//
//     vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
//     reader->SetFileName(file_name.c_str());
//     reader->Update();
//     vtkUnstructuredGrid *vtk_grid = reader->GetOutput();
//     vtk_grid->Update();
//
//     for (vtkIdType i = 0; i < vtk_grid->GetNumberOfPoints(); ++i)
//     {
//         double * point = vtk_grid->GetPoint(i);
//         vector3_type v(point[0], point[1], point[2]);
//         mesh.add_vertex(v);
//     }
//
//     for (vtkIdType i = 0; i < vtk_grid->GetNumberOfCells(); ++i)
//     {
//         vtkCell* cell = vtk_grid->GetCell(i);
//         int size = cell->GetNumberOfPoints();
//         vtkIdType *points = cell->GetPointIds()->GetPointer(0);
//         std::list<trimesh::TriMesh<>::vertex_handle> handles;
//
//         for (int j = 0; j < size; ++j)
//             handles.push_back(mesh.get_vertex_handle(points[j]));
//
//         mesh.add_face(handles.begin(), handles.end());
//     }
//
//     reader->Delete();
// }

    template<typename vector_type>
    bool read_vtk ( int grid_form, const std::string &file_name, vector_type &p, vector_type &v, vector_type &f )
    {
        typedef typename vector_type::value_type real_type;
        typedef typename vector_type::size_type  index_type;

        switch ( grid_form )
        {

            case 0:
            {
                vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
                reader->SetFileName ( file_name.c_str() );
                reader->Update();
                vtkUnstructuredGrid *vtk_grid = reader->GetOutput();
                vtk_grid->Update();

                vtkDoubleArray* velocities = static_cast<vtkDoubleArray*> ( vtk_grid->GetPointData()->GetVectors ( "velocities" ) );
                vtkDoubleArray* forces = static_cast<vtkDoubleArray*> ( vtk_grid->GetPointData()->GetVectors ( "forces" ) );

                index_type size = 3*vtk_grid->GetNumberOfPoints();
                p.resize ( size );
                v.resize ( size );
                f.resize ( size );
                for ( vtkIdType i = 0,j = 0; i < vtk_grid->GetNumberOfPoints(), j < size; ++i, j+=3 )
                {
                    // Set Positions
                    real_type * tuple = vtk_grid->GetPoint ( i );
                    p ( j ) = tuple[0];
                    p ( j+1 ) = tuple[1];
                    p ( j+2 ) = tuple[2];
                    // Set Velocities
                    velocities->GetTupleValue ( i, tuple );
                    v ( j ) = tuple[0];
                    v ( j+1 ) = tuple[1];
                    v ( j+2 ) = tuple[2];
                    // Set Forces
                    forces->GetTupleValue ( i, tuple );
                    f ( j ) = tuple[0];
                    f ( j+1 ) = tuple[1];
                    f ( j+2 ) = tuple[2];
                }
                reader->Delete();

                break;
            }

            case 1:
            {
                vtkXMLStructuredGridReader *reader = vtkXMLStructuredGridReader::New();
                reader->SetFileName ( file_name.c_str() );
                reader->Update();
                vtkStructuredGrid *vtk_grid = reader->GetOutput();
                vtk_grid->Update();

                vtkDoubleArray* velocities = static_cast<vtkDoubleArray*> ( vtk_grid->GetPointData()->GetVectors ( "velocities" ) );
                vtkDoubleArray* forces = static_cast<vtkDoubleArray*> ( vtk_grid->GetPointData()->GetVectors ( "forces" ) );

                index_type size = 3*vtk_grid->GetNumberOfPoints();
                p.resize ( size );
                v.resize ( size );
                f.resize ( size );
                for ( vtkIdType i = 0,j = 0; i < vtk_grid->GetNumberOfPoints(), j < size; ++i, j+=3 )
                {
                    // Set Positions
                    real_type * tuple = vtk_grid->GetPoint ( i );
                    p ( j ) = tuple[0];
                    p ( j+1 ) = tuple[1];
                    p ( j+2 ) = tuple[2];
                    // Set Velocities
                    velocities->GetTupleValue ( i, tuple );
                    v ( j ) = tuple[0];
                    v ( j+1 ) = tuple[1];
                    v ( j+2 ) = tuple[2];
                    // Set Forces
                    forces->GetTupleValue ( i, tuple );
                    f ( j ) = tuple[0];
                    f ( j+1 ) = tuple[1];
                    f ( j+2 ) = tuple[2];
                }
                reader->Delete();

                break;
            }

            default:
            {
                std::cerr << "Grid type not implemented." << std::endl;
                return 0;
            }
        }

        return 1;

    }



}
#endif
// kate: indent-mode cstyle; space-indent on; indent-width 4;
