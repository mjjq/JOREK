Name = 'JorekReadH5'
Label = 'Jorek Read H5 file'
Help = 'Read a jorek restart file in H5 format'
Extension = 'h5'
FileDescription = 'JOREK hdf5 files'

NumberOfInputs = 0
OutputDataType = 'vtkUnstructuredGrid'
Properties = dict(
    Variables = "1,2,3,4,5,6"
)

#from paraview import vtk
import numpy as np
from vtk import vtkPointData, vtkDataArray, vtkUnstructuredGrid, vtkPoints, \
     vtkIdList
import os
import jorek_read_h5 as jorek
import vtk.numpy_interface as vtknp

def RequestData():
    pdo = self.GetOutput()
    #output = self.getUnstructuredGridOutput()

    # Read the h5 file
    fields = jorek.fields(Variables)
    fields.read(FileName)

    # Create 2D interpolation
    (x, values) = fields.interp_to_mesh()

    numPoints = x.shape[0]*x.shape[1]*x.shape[2]
    RZ = np.reshape(x, (numPoints, 2))
    #output.PointData['coordsX'] = RZ[:,1]
    #output.PointData['coordsY'] = RZ[:,2]
    pts = vtk.vtkPoints()
    for i in range(numPoints):
        pts.InsertNextPoint(RZ[i,0], RZ[i,1], 0)
    pdo.SetPoints(pts)
