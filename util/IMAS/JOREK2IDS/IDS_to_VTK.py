import numpy as np
import getpass
import argparse
import vtk

from vtk.util import numpy_support as npvtk
from idsUtilities import basicIDS, readIDS

prec=np.float32
vtk_prec=vtk.VTK_FLOAT


parser = argparse.ArgumentParser(description="Convert IMAS MHD IDS to VTK file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--shot", type=int, default=1, help="Shot number")
parser.add_argument("-r", "--run", type=int, default=1, help="Run number")
parser.add_argument("-u", "--user", type=str, default=getpass.getuser(),
                    help="Location of ~$USER/public/imasdb")
parser.add_argument("-d", "--database", type=str, default="smiter", help="Database name under public/imasdb/")
parser.add_argument("-o", "--occurrence", type=int, default=0, help="Occurrence number")
parser.add_argument("vtkfile", metavar='jorek.vtk', nargs='?', help="Resulting VTK filename", default="jorek.vtk")
args = parser.parse_args()

b_ids = basicIDS(args.shot, args.run, args.user, args.database)
r_ids = readIDS(args.shot, args.run, args.user, args.database)
r_ids.getGGD("mhd")


#grid and conectivity
xyz0 = r_ids.grid_ggd.array[0].space.array[0].objects_per_dimension.array[0].object.array
xyz = np.empty((len(xyz0), 3))
for i in range(len(xyz0)):
    xyz[i, :] = xyz0[i].geometry

"""ien = np.empty((len(xyz0) // 4, 4))
for j in range(len(xyz0) // 4):
    s = int(j * 4)
    ien[j, :] = np.array([s, s + 2, s + 3, s + 1])"""

ien0 = np.array(r_ids.grid_ggd.array[0].space.array[0].objects_per_dimension.array[2].object.array)
ien = np.empty((np.shape(ien0)[0], np.shape(ien0[0].nodes)[0]))
for i in range(np.shape(ien0)[0]):
    ien[i, :] = np.array(ien0[i].nodes) - 1
ien = np.insert(ien, 0, np.shape(ien0[0].nodes)[0], axis=1)

#values
v =[]
if r_ids.ggd.array[0].electrons.temperature.array != []:
    v = np.array(r_ids.ggd.array[0].electrons.temperature.array[0].values)

#write vtk
etype = vtk.VTK_QUAD
output = vtk.vtkUnstructuredGrid()

pcoords = npvtk.numpy_to_vtk(xyz, deep=True, array_type=vtk_prec)
p = vtk.vtkPoints()
p.SetData(pcoords)
output.SetPoints(p)

c = vtk.vtkCellArray()
c.SetCells(ien.shape[0], npvtk.numpy_to_vtk(ien, deep=True, array_type=vtk.VTK_ID_TYPE))
output.SetCells(etype, c)

if v != []:
    tmp = npvtk.numpy_to_vtk(v, deep=True, array_type=vtk_prec)
    tmp.SetName("T")
    output.GetPointData().AddArray(tmp)

writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileName(args.vtkfile)
writer.SetInputData(output)
writer.Write()
exit(0)
