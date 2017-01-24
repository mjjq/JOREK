Name = 'JorekReadH5'
Label = 'Jorek Read H5 file'
Help = 'Read a jorek restart file in H5 format'
Extension = 'h5'
FileDescription = 'JOREK hdf5 files'

NumberOfInputs = 0
OutputDataType = 'vtkUnstructuredGrid'
Properties = dict(
    Variables = "1,2,3,4,5,6",
    Number_of_planes = 16,
    Number_of_subdivisions = 5,
    phi_range_in_pi = [0.0, 0.5],
)


# from paraview import vtk # is done automatically
def RequestData(self):
    import jorek_read_h5 as jorek
    output = self.GetUnstructuredGridOutput()

    # Read the h5 file
    #fields = jorek.fields()
    #fields.read(FileName, variables=Variables)
    
    #pdo = fields.to_vtk(phi=0)
    #return pdo
