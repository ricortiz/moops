/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDiskSource.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSurfaceSource.h"

#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

vtkStandardNewMacro(vtkSurfaceSource);

vtkSurfaceSource::vtkSurfaceSource() 
{ 
  this->SetNumberOfInputPorts(0);   
};

vtkSurfaceSource::vtkSurfaceSource(int i)
{
    this->SetNumberOfInputPorts(i);
};

void vtkSurfaceSource::SetPositions(vtkPoints *points)
{
    m_points = points;
}

void vtkSurfaceSource::SetCells(vtkCellArray *cells)
{
    m_cells = cells;
}

int vtkSurfaceSource::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
{
    // get the info object
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the ouptut
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Update ourselves and release memory
    //
    if(!output)
    {
        output->SetPoints(m_points);
        output->SetPolys(m_cells);
    }
    return 1;
}

void vtkSurfaceSource::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}
