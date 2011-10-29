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

// Label for the type of grid you want to read
enum {unstructured, structured};

/// @brief This function reads VTK unstructured(.vtu) and structured(.vts) grids
/// @param grid_form Can be either unstructured or structured
/// @param file_name Should be a valid VTK xml format file (either vtu or vts)
/// @param p upon return holds position coordinates of all points
/// @param v upon return holds velocities of all points
/// @param f upon return holds forces of all points
template<typename vector_type>
bool read_vtk ( int grid_form, const std::string &file_name, vector_type &p, vector_type &v, vector_type &f )
{
    typedef typename vector_type::value_type real_type;
    typedef typename vector_type::size_type  index_type;

    switch ( grid_form )
    {

    case unstructured:
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

    case structured:
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
