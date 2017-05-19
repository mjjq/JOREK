Name = 'JorekReadH5Proj'
Label = 'Jorek Read Projection H5 file'
Help = 'Read a jorek projection file in H5 format'
Extension = 'h5'
FileDescription = 'JOREK Projection hdf5 files'

import h5py
import inspect
from os.path import dirname
PythonPaths = [dirname(dirname(inspect.getfile(h5py)))]

NumberOfInputs = 0
OutputDataType = 'vtkUnstructuredGrid'
Properties = dict(
    Variables = "1",
    Number_of_planes = 1,
    Number_of_subdivisions = 4,
    phi_range_in_pi = [0.0, 1.0],
    Quadratic = False,
    Exclude_n0_mode = False,
    central_mass = 2.0, # Convert times to SI units using this norm
    central_density = 1.0 # [1d-20 m^-3]
)


# from paraview import vtk # is done automatically
def RequestData(self):
    #include 'proj_read_h5.py' # Don't delete: include content of this file here

    def GetUpdateTimestep(algorithm):
        """Returns the requested time value, or None if not present"""
        executive = algorithm.GetExecutive()
        outInfo = executive.GetOutputInformation(0)
        return outInfo.Get(executive.UPDATE_TIME_STEP()) \
                if outInfo.Has(executive.UPDATE_TIME_STEP()) else None
    # Get the current timestep
    req_time = GetUpdateTimestep(self)

    # Read the timestep info from the files
    xtime = [h5py.File(fname).get("t_now")[0] for fname in FileNames]

    # 4 possibilities here:
    # After last step: return last file
    # Before first step: return first file
    # Very close match: return that file
    # Between 2 files: interpolate
    interp = False
    try:
        index = np.isclose(xtime, req_time).tolist().index(True)
        # We have a match, read and return data for this filename
    except ValueError:
        # Check for other 3 cases
        if (np.count_nonzero(req_time < xtime) == 0):
            # After last file
            index = len(FileNames)-1
        elif (np.count_nonzero(req_time > xtime) == 0):
            # Before first file
            index = 0
        else:
            interp = True
            index = (req_time > xtime).tolist().index(True)
            f = (req_time - xtime[index-1])/(xtime[index] - xtime[index-1])
            # how much of second to take == 1-how much of first to take

    # Read the h5 file
    if (not hasattr(self, 'f')):
        self.f = fields()
    if (interp):
        self.f.read(FileNames[index], variables=Variables, file_prev=FileNames[index-1],
               interp_fraction=f)
    else:
        self.f.read(FileNames[index], variables=Variables)

    output = self.f.to_vtk(n_sub=Number_of_subdivisions, phi=phi_range_in_pi,
                      n_plane=Number_of_planes, without_n0_mode=Exclude_n0_mode,
                      output=self.GetUnstructuredGridOutput(), quadratic=Quadratic)
    return output

"""
See paraview guide 13.2.2
"""
def RequestInformation(self):
    import numpy as np
    import h5py
    def setOutputTimesteps(algorithm):
        executive = algorithm.GetExecutive()
        outInfo = executive.GetOutputInformation(0)

        # Read the timestep info from the files
        xtime = [h5py.File(fname).get("t_now")[0] for fname in FileNames]

        outInfo.Remove(executive.TIME_STEPS())
        for i in range(len(FileNames)):
            outInfo.Append(executive.TIME_STEPS(), xtime[i])

        # Remove time range info
        outInfo.Remove(executive.TIME_RANGE())
        outInfo.Append(executive.TIME_RANGE(), xtime[0])
        outInfo.Append(executive.TIME_RANGE(), xtime[-1])

    setOutputTimesteps(self)
