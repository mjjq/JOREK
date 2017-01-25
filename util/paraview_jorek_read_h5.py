Name = 'JorekReadH5'
Label = 'Jorek Read H5 file'
Help = 'Read a jorek restart file in H5 format'
Extension = 'h5'
FileDescription = 'JOREK hdf5 files'

import h5py
import inspect
from os.path import dirname
PythonPaths = [dirname(dirname(inspect.getfile(h5py)))]

NumberOfInputs = 0
OutputDataType = 'vtkUnstructuredGrid'
Properties = dict(
    Variables = "1,2,3,4,5,6",
    Number_of_planes = 1,
    Number_of_subdivisions = 5,
    phi_range_in_pi = [0.0, 0.0],
    Exclude_n0_mode = False,
)


# from paraview import vtk # is done automatically
def RequestData(self):
    #include 'jorek_read_h5.py' # Don't delete: include content of this file here

    def GetUpdateTimestep(algorithm):
        """Returns the requested time value, or None if not present"""
        executive = algorithm.GetExecutive()
        outInfo = executive.GetOutputInformation(0)
        return outInfo.Get(executive.UPDATE_TIME_STEP()) \
                if outInfo.Has(executive.UPDATE_TIME_STEP()) else None
    # Get the current timestep
    req_time = GetUpdateTimestep(self)

    req_time = round(req_time)
    if (req_time < 0):
        req_time = 0
    elif (req_time >= len(FileNames)):
        req_time = len(FileNames)-1

    fname = FileNames[req_time]
    print(fname)
    # Read the h5 file
    f = fields()
    f.read(fname, variables=Variables)
    
    output = fields.to_vtk(n_sub=Number_of_subdivisions, phi=phi_range_in_pi,
                           n_plane=Number_of_planes, without_n0_mode=Exclude_n0_mode,
                           output=self.GetUnstructuredGridOutput())
    return output

"""
See paraview guide 13.2.2
"""
def RequestInformation(self):
    def setOutputTimesteps(algorithm):
        executive = algorithm.GetExecutive()
        outInfo = executive.GetOutputInformation(0)

        outInfo.Remove(executive.TIME_STEPS())
        for timestep in range(len(FileNames)):
            outInfo.Append(executive.TIME_STEPS(), timestep)

        # Remove time range info
        outInfo.Remove(executive.TIME_RANGE())
        outInfo.Append(executive.TIME_RANGE(), 0)
        outInfo.Append(executive.TIME_RANGE(), len(FileNames))

    setOutputTimesteps(self)
