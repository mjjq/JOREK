import vtk
import numpy as np
import jorek_read_h5 as jorek

f = jorek.fields()
f.read('/tmp/jorek_restart.h5', variables='1')
grid = f.to_vtk()

 
# Visualize
mapper = vtk.vtkDataSetMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(grid)
else:
    mapper.SetInputData(grid)

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

#WritePNG(renderWindowInteractor.GetRenderWindow(), "Cell3DDemonstration.png")
