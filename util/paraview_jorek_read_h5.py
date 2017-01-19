Name = 'JorekReadH5'
Label = 'Jorek Read H5 file'
Help = 'Read a jorek restart file in H5 format'

NumberOfInputs = 0
OutputDataType = 'vtkUnstructuredGrid'
Properties = dict(
    Variables = "1,2,3,4,5,6",
    FileName = "jorek_restart.h5"
)

#from paraview import vtk
import numpy as np
from vtk import vtkPointData, vtkDataArray, vtkUnstructuredGrid, vtkPoints, \
     vtkIdList
import os
import jorek_read_h5 as jorek
import vtk.numpy_interface as vtknp

def RequestData():
    output = self.getUnstructuredGridOutput()

    # Read the h5 file
    fields = jorek.fields(Variables)
    fields.read(FileName)

    # Create 2D interpolation
    (x, values) = fields.interp_to_mesh()

    numPoints = x.shape[0]*x.shape[1]*x.shape[2]
    RZ = np.reshape(x, (numPoints, 2))
    output.PointData['coordsX'] = RZ[:,1]
    output.PointData['coordsY'] = RZ[:,2]
