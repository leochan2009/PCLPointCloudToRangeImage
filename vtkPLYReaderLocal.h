/*=========================================================================
 
 Program:   Visualization Toolkit
 Module:    vtkPLYReaderLocal.h
 
 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.
 
 =========================================================================*/
// .NAME vtkPLYReaderLocal - read Stanford University PLY polygonal file format
// .SECTION Description
// vtkPLYReaderLocal is a source object that reads polygonal data in
// Stanford University PLY file format (see
// http://graphics.stanford.edu/data/3Dscanrep). It requires that
// the elements "vertex" and "face" are defined. The "vertex" element
// must have the properties "x", "y", and "z". The "face" element must
// have the property "vertex_indices" defined. Optionally, if the "face"
// element has the properties "intensity" and/or the triplet "red",
// "green", and "blue"; these are read and added as scalars to the
// output data.

// .SECTION See Also
// vtkPLYWriter

#ifndef vtkPLYReaderLocal_h
#define vtkPLYReaderLocal_h

#include "vtkAbstractPolyDataReader.h"

class vtkPLYReaderLocal : public vtkAbstractPolyDataReader
{
public:
  vtkTypeMacro(vtkPLYReaderLocal,vtkAbstractPolyDataReader);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  // Description:
  // Construct object with merging set to true.
  static vtkPLYReaderLocal *New();
  
  // Description:
  // A simple, non-exhaustive check to see if a file is a valid ply file.
  static int CanReadFile(const char *filename);
  vtkPolyData * GetPolyData(){return output;};

protected:
  vtkPLYReaderLocal();
  ~vtkPLYReaderLocal();
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
private:
  vtkPLYReaderLocal(const vtkPLYReaderLocal&);  // Not implemented.
  void operator=(const vtkPLYReaderLocal&);  // Not implemented.
  vtkPolyData *output;
};

#endif
