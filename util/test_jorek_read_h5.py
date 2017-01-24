from __future__ import print_function
import vtk
import numpy as np
import jorek_read_h5 as jorek
import sys

f = jorek.fields()
f.read('/tmp/jorek_restart.h5', variables='1')
grid = f.to_vtk(phi=0)
#grid = f.to_vtk(phi=[0,np.pi/2])

#sys.exit(0)


# Visualize
mapper = vtk.vtkDataSetMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(grid)
else:
    mapper.SetInputData(grid)

# Create color map based on range of first var
data_range = grid.GetPointData().GetAbstractArray(0).GetRange()
colormap = vtk.vtkLookupTable()
colormap.SetHueRange(data_range[0], data_range[1])
colormap.Build()
mapper.SetColorModeToDefault()
mapper.SetScalarRange(data_range)
mapper.SetScalarVisibility(True)
mapper.SetLookupTable(colormap)

actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetPointSize(20)

renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

renderer.AddActor(actor)

renderWindow.Render()
renderWindowInteractor.Start()
