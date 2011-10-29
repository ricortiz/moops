/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSurfaceSource.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSurfaceSource - create a disk with hole in center
// .SECTION Description
// vtkSurfaceSource creates a polygonal disk with a hole in the center. The 
// disk has zero height. The user can specify the inner and outer radius
// of the disk, and the radial and circumferential resolution of the 
// polygonal representation. 
// .SECTION See Also
// vtkLinearExtrusionFilter

#ifndef __vtkSurfaceSource_h
#define __vtkSurfaceSource_h

class vtkCellArray;
class vtkPoints;

#include <vtkPolyDataAlgorithm.h>

class VTK_GRAPHICS_EXPORT vtkSurfaceSource : public vtkPolyDataAlgorithm 
{
public:
  static vtkSurfaceSource *New();
  vtkTypeMacro(vtkSurfaceSource,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  void SetPositions(vtkPoints *surface);  
  void SetCells(vtkCellArray *cells);

protected:
    vtkSurfaceSource();
    vtkSurfaceSource(int);
  ~vtkSurfaceSource() {};

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
  vtkPoints    *m_points;
  vtkCellArray *m_cells;

private:
  vtkSurfaceSource(const vtkSurfaceSource&);  // Not implemented.
  void operator=(const vtkSurfaceSource&);  // Not implemented.
};

#endif
