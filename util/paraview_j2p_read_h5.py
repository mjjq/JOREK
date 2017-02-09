Name = 'J2PReadH5'
Label = 'Jorek Particles Read H5 file'
Help = 'Read a jorek particle restart file in H5 format'
Extension = 'h5'
FileDescription = 'JOREK particle hdf5 files'

import h5py
import inspect
from os.path import dirname
PythonPaths = [dirname(dirname(inspect.getfile(h5py)))]

NumberOfInputs = 0
OutputDataType = 'vtkPolyData'
Properties = dict(
    group='1',
    toroidal=True,
    select='', # comma-separated list of particle numbers, 0-indexed
)


# from paraview import vtk # is done automatically
def RequestData(self):
    import h5py
    import numpy as np
    import vtk
    from vtk.util import numpy_support as npvtk
    import re
    time_re = r"\d+\.?\d*"

    def GetUpdateTimestep(algorithm):
        """Returns the requested time value, or None if not present"""
        executive = algorithm.GetExecutive()
        outInfo = executive.GetOutputInformation(0)
        return outInfo.Get(executive.UPDATE_TIME_STEP()) \
                if outInfo.Has(executive.UPDATE_TIME_STEP()) else None
    # Get the current timestep
    req_time = GetUpdateTimestep(self)

    # Find the closest two files
    times = np.asarray([float(re.search(time_re, f).group()) for f in FileNames])
    # 4 possibilities here:
    # After last step: return last file
    # Before first step: return first file
    # Very close match: return that file
    # Between 2 files: interpolate
    interp = False
    try:
        index = np.isclose(times, req_time).tolist().index(True)
        # We have a match, read and return data for this filename
    except ValueError:
        # Check for other 3 cases
        if (np.count_nonzero(req_time < times) == 0):
            # After last file
            index = len(FileNames)-1
        elif (np.count_nonzero(req_time > times) == 0):
            # Before first file
            index = 0
        else:
            interp = True
            index = (req_time > times).tolist().index(True)
            f = (req_time - times[index-1])/(times[index] - times[index-1])
            # how much of second to take == 1-how much of first to take

    gid = int(group)
    if (len(select) > 0):
        sel = np.s_[[int(s) for s in select.split(',')]]
    else:
        sel = np.s_[:]

    # Read the h5 file(s)
    with h5py.File(FileNames[index]) as hf:
        x      = hf.get('groups/%03d/x'%gid)[sel,:]
        weight = hf.get('groups/%03d/weight'%gid)[sel]
        q      = hf.get('groups/%03d/q'%gid)[sel]
    if (interp):
        with h5py.File(FileNames[index-1]) as hf:
            x      = f*x      + (1.0-f)*hf.get('groups/%03d/x'%gid)[sel,:]
            weight = f*weight + (1.0-f)*hf.get('groups/%03d/weight'%gid)[sel,:]
            q      = f*q      + (1.0-f)*hf.get('groups/%03d/q'%gid)[sel,:]

    if (not toroidal):
        pcoords = npvtk.numpy_to_vtk(x, deep=True, array_type=vtk.VTK_FLOAT)
    else:
        tmp = np.stack((x[:,0]*np.cos(x[:,2]),
                        x[:,1],
                        x[:,0]*np.sin(x[:,2])), axis=-1)
        pcoords = npvtk.numpy_to_vtk(tmp, deep=True, array_type=vtk.VTK_FLOAT)

    points = vtk.vtkPoints()
    points.SetData(pcoords)
    output.SetPoints(points)

    val = npvtk.numpy_to_vtk(weight, deep=True, array_type=vtk.VTK_FLOAT)
    val.SetName("weight")
    output.GetPointData().AddArray(val)

    val = npvtk.numpy_to_vtk(q, deep=True, array_type=vtk.VTK_FLOAT)
    val.SetName("q")
    output.GetPointData().AddArray(val)

    return output

"""
See paraview guide 13.2.2
"""
def RequestInformation(self):
    import re
    time_re = r"\d+\.?\d*"
    def setOutputTimesteps(algorithm):
        executive = algorithm.GetExecutive()
        outInfo = executive.GetOutputInformation(0)

        outInfo.Remove(executive.TIME_STEPS())
        for filename in FileNames:
            # keep only the numbers and dots and remove the last dot
            timestep = float(re.search(time_re, filename).group())
            outInfo.Append(executive.TIME_STEPS(), timestep)

        # Remove time range info
        outInfo.Remove(executive.TIME_RANGE())
        outInfo.Append(executive.TIME_RANGE(), float(re.search(time_re, FileNames[0]).group()))
        outInfo.Append(executive.TIME_RANGE(), float(re.search(time_re, FileNames[-1]).group()))

    setOutputTimesteps(self)
