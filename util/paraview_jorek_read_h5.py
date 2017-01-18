Name = 'JorekReadH5'
Label = 'Jorek Read H5 file'
Help = 'Read a jorek restart file in H5 format'

NumberOfInputs = 0
OutputDataType = 'vtkUnstructuredGrid'
properties = dict(
    Variables = "1,2,3,4,5,6",
    FileName = "jorek_restart.h5"
)

from paraview import vtk
import os
import jorek_read_h5 as jorek

output = self.getUnstructuredGridOutput()

# Read the h5 file
fields = jorek.fields(Variables)
fields.read(FileName)

# 
