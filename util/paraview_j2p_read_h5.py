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

    def GetUpdateTimestep(algorithm):
        """Returns the requested time value, or None if not present"""
        executive = algorithm.GetExecutive()
        outInfo = executive.GetOutputInformation(0)
        return outInfo.Get(executive.UPDATE_TIME_STEP()) \
                if outInfo.Has(executive.UPDATE_TIME_STEP()) else None
    # Get the current timestep
    req_time = GetUpdateTimestep(self)

    req_time = int(round(req_time))
    if (req_time < 0):
        req_time = 0
    elif (req_time >= len(FileNames)):
        req_time = len(FileNames)-1

    fname = FileNames[req_time]

    gid = int(group)
    if (len(select) > 0):
        sel = np.s_[[int(s) for s in select.split(',')]]
    else:
        sel = np.s_[:]

    # Read the h5 file
    with h5py.File(fname) as hf:
        x      = hf.get('groups/%03d/x'%gid)
        weight = hf.get('groups/%03d/weight'%gid)
        q      = hf.get('groups/%03d/q'%gid)
        i_elm  = hf.get('groups/%03d/i_elm'%gid)

        if (not toroidal):
            pcoords = npvtk.numpy_to_vtk(x[sel,:], deep=True, array_type=vtk.VTK_FLOAT)
        else:
            tmp = np.stack((x[sel,0]*np.cos(x[sel,2]),
                            x[sel,1],
                            x[sel,0]*np.sin(x[sel,2])), axis=-1)
            pcoords = npvtk.numpy_to_vtk(tmp, deep=True, array_type=vtk.VTK_FLOAT)
        points = vtk.vtkPoints()
        points.SetData(pcoords)
        output.SetPoints(points)

        val = npvtk.numpy_to_vtk(weight[sel], deep=True, array_type=vtk.VTK_FLOAT)
        val.SetName("weight")
        output.GetPointData().AddArray(val)

        val = npvtk.numpy_to_vtk(q[sel], deep=True, array_type=vtk.VTK_FLOAT)
        val.SetName("q")
        output.GetPointData().AddArray(val)

        val = npvtk.numpy_to_vtk(i_elm[sel], deep=True, array_type=vtk.VTK_INT)
        val.SetName("i_elm")
        output.GetPointData().AddArray(val)

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
