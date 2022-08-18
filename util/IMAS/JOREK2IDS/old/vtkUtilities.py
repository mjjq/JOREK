#!/usr/bin/env python

#  Name : jorekHDF5toIDS.py
#
#  Description :
#           A collection of routines for handling VTK files.
#
#           Requirements:
#           - pip3 install --user vtk
#           - pip3 install --user pyvtk
#           - pip3 install --user pyqt5
#
#  Author :
#         Dejan Penko
#  E-mail :
#         dejan.penko@lecad.fs.uni-lj.si
#
#****************************************************
#     Copyright(c) 2020- D. Penko

import vtk
from vtk.util import numpy_support as VN
import logging
import sys
import numpy as np

from PyQt5 import QtWidgets

class convert():

    def __init__(self):
        pass

    def binary2ascii(self, filePath):
        reader = vtk.vtkDataSetReader()
        # reader = vtk.vtkUnstructuredGridReader()
        # reader = vtk.vtkPolyDataReader()
        reader.SetFileName(filePath)
        reader.Update()
        reader.UpdateDataObject()
        readerData  = reader.GetOutput()

        print(readerData)

        w = vtk.vtkDataSetWriter()
        # w.SetPrecision(10)
        # w.SetScientific()
        w.SetInputData( readerData )
        w.SetFileName('foo-ascii.vtk')
        w.Write()

class vtkReadUtilities:
    def __init__(self, fileType='ASCII'):

        # NOTE: Python vtk library has troubles reading all field data from
        # binary files!

        # Empty reader data object
        self.readerData = None
        # Get path to file with PyQt dialog
        self.filePath = self.getVtkFileDialog()
        # Set readerData object
        self.setReaderData(self.filePath, fileType)

        self.readerData.points_geo = []
        self.readerData.obj_0D_list = []
        self.readerData.obj_1D_list = []
        self.readerData.obj_2D_list = []
        self.readerData.obj_3D_list = []

        self.readerData.obj_0D_type = None
        self.readerData.obj_1D_type = None
        self.readerData.obj_2D_type = None
        self.readerData.obj_3D_type = None

        self.readerData.field_dict = {}

    def getVtkFileDialog(self):
        """Run as standalone application.
        """
        # Set QApplication
        app = QtWidgets.QApplication([])
        # Set options
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        # Set file dialog
        # Note: arguments are:  - parent (None),
        #                       - window title
        #                       - default file name
        #                       - file types selectable options
        #                       - options
        filePath, _ = QtWidgets.QFileDialog \
            .getOpenFileName(None,
                             "Select VTK file",
                             "",
                             "All Files (*);;VTK Files (*.vtk)",
                             options=options)
        # Close application
        app.exit()

        if filePath == '':
            logging.error("No VTK file was selected! Exiting the program.")
            sys.exit()

        print("Selected VTK file (full path): ", filePath)

        return filePath

    def setReaderData(self, filePath, fileType="ASCII"):
        """Set readerData object.

        Arguments:
            filePath (str) : Full path to file.
            fileTyoe (str) : Type of file, ASCII (default) or XML.
        """

        print("Setting VTK reader")

        if fileType == "ASCII":
            self.reader = vtk.vtkUnstructuredGridReader()    # For ASCII/BINARY
            # self.reader = vtk.vtkDataSetReader()    # For ASCII/BINARY
        elif fileType == "XML":
            self.reader = vtk.vtkXMLUnstructuredGridReader() # For XML files
        else:
            logging.warning(" VTK file type not properly set.")
        # Notes:
        # Get the type of file (ASCII or BINARY). 1 --> ASCII
        #                                         2 --> BINARY
        fileTypeID = self.reader.GetFileType()

        if fileTypeID == 1:
            print("FileType: ASCII")
        elif fileTypeID == 2:
            print("FileType: BINARY")
        else:
            print("FileType: INVALID!")

        # Read file and set reader data object
        self.reader.SetFileName(filePath)
        self.reader.Update()
        # Set output
        self.readerData = self.reader.GetOutput()

    def getReaderData(self):
        return self.readerData

    def readVtkUnstructuredGrid(self):
        """ Extract points geometry/coordinates and connectivity arrays for 0D,
        1D, 2D and 3D elements. Works both for ASCII and BINARY file type.

        arguments:
            filePath (str) : Full path (including file name) to ASCII/BINARY
                             VTK file.
        """

        print("START: Reading VTK grid")

        # Empty lists for point geometry/coordinates
        points_geo = []
        # Empty lists for n-D objects/cells
        obj_0D_list = []    # vertices
        obj_1D_list = []    # edges
        obj_2D_list = []    # 2D cells (polys)
        obj_3D_list = []    # 3D volumes

        # Read points geometry in 3D space
        # POINTS
        points = self.readerData.GetPoints()
        logging.info("First point: ", self.readerData.GetPoint(0))
        array_points = points.GetData()
        points_geo = VN.vtk_to_numpy(array_points)

        # Note: the precision is fine. only numpy presentation shown only 5
        #       decimal numbers
        print("Nodes: ", points_geo)

        print("Number of nodes: ", len(points_geo))
        print("*************************************")
        # POINTS =/= VERTICES (NODES)!!! VERTICES ARE SET AS TYPE OF CELLS


        # num_0D_obj = self.readerData.GetNumberOfVerts()
        # print("num_0D_obj: ", num_0D_obj)
        # num_1D_obj = self.readerData.GetNumberOfLines()
        # print("num_1D_obj: ", num_1D_obj)
        # num_2D_obj = self.readerData.GetNumberOfPolys()
        # print("num_2D_obj: ", num_2D_obj)
        num_all_obj = self.readerData.GetNumberOfCells()
        print("num_all_obj: ", num_all_obj)
        max_cell_size = self.readerData.GetMaxCellSize()
        print("max_cell_size: ", max_cell_size)

        for i in range(0, num_all_obj):
            obj = self.readerData.GetCell(i)   # Traverses through all cells: vertices,
                                  # then lines, then polys
            obj_type = self.readerData.GetCellType(i)  # VERTIEX:  type ID = 1
                                                       # LINE:     type ID = 3
                                                       # POLY:     type ID = 5

            """
            VTK descriptor| VTK name                  | Finite Element type
            --------------|---------------------------|-------------------------
            1               VTK_VERTEX                  Vertex
            2               VTK_POLY_VERTEX             Vertex
            3               VTK_LINE                    Edge Lagrange P1
            5               VTK_TRIANGLE                Triangle Lagrange P1
            8               VTK_PIXEL                   Quadrilateral Lagrange P1
            9               VTK_QUAD                    Quadrilateral Lagrange P1
            10              VTK_TETRA                   Tetrahedron Lagrange P1
            11              VTK_VOXEL                   Hexahedron Lagrange P1
            12              VTK_HEXAHEDRON              Hexahedron Lagrange P1
            13              VTK_WEDGE                   Wedge Lagrange P1
            21              VTK_QUADRATIC_EDGE          Edge Lagrange P2
            22              VTK_QUADRATIC_TRIANGLE      Triangle Lagrange P2
            23              VTK_QUADRATIC_QUAD          Quadrilateral Lagrange P2
            24              VTK_QUADRATIC_TETRA         Tetrahedron Lagrange P2
            25              VTK_QUADRATIC_HEXAHEDRON    Hexahedron Lagrange P2
            """
            if obj_type == vtk.VTK_VERTEX: # 1
                obj_0D_list.append(obj.GetPointId(0))
            elif obj_type == vtk.VTK_LINE: # 3
                obj_1D_list.append((obj.GetPointId(0), obj.GetPointId(1)))
            elif obj_type == vtk.VTK_TRIANGLE: # 5
                obj_2D_list.append((obj.GetPointId(0), obj.GetPointId(1),
                                    obj.GetPointId(2)))
            elif obj_type == vtk.VTK_QUAD: # 9
                obj_2D_list.append((obj.GetPointId(0), obj.GetPointId(1),
                                    obj.GetPointId(2), obj.GetPointId(3)))
            elif obj_type == vtk.VTK_TETRA: # 10
                obj_3D_list.append((obj.GetPointId(0), obj.GetPointId(1),
                                    obj.GetPointId(2), obj.GetPointId(3)))

        # Print summary
        print("Number of points: ", len(points_geo))
        print("Number of 1D elements: ", len(obj_1D_list))
        print("Number of 2D elements: ", len(obj_2D_list))
        print("Number of 3D elements: ", len(obj_3D_list))
        if len(points_geo) > 0:
            print("First point: ", points_geo[0])
        else:
            print("First point: None")
        if len(obj_0D_list) > 0:
            print("First obj_0D: ", obj_0D_list[0])
        else:
            print("First obj_0D: None")
        if len(obj_1D_list) > 0:
            print("First obj_1D: ", obj_1D_list[0])
        else:
            print("First obj_1D: None")
        if len(obj_2D_list) > 0:
            print("First obj_2D: ", obj_2D_list[0])
        else:
            print("First obj_2D: None")
        if len(obj_3D_list) > 0:
            print("First obj_3D: ", obj_3D_list[0])
        else:
            print("First obj_3D: None")

        print("vtk_read_unstructuredgrid FINISHED")

        # Set values to readerData
        self.readerData.points_geo = np.array(points_geo)
        self.readerData.obj_0D_list = np.array(obj_0D_list)
        self.readerData.obj_1D_list = np.array(obj_1D_list)
        self.readerData.obj_2D_list = np.array(obj_2D_list)
        self.readerData.obj_3D_list = np.array(obj_3D_list)

    def readVtkDataFields(self, fileType="ASCII"):
        """ Extracts all field data arrays.
        Note: only POINT DATA are currently extracted!
        TODO: read CELL DATA.
        Populates rvfd object "field_dict" dictionary. Its definition:

            field_dict = {"array_label_i"  : [i_1, i_2, ..., i_n],
                 "array_label_j"  : [j_1, j_2, ..., j_n],
                 ...
                 "array_label_m"  : [m_1, m_2, ..., m_n]}

        arguments:
            filePath (str) : Full path (including file name) to ASCII TK file.
            fileType (str) : File type (available options: "ASCII" and "XML")
        """

        print("START: Reading VTK fields")

        # FIELDS RELATED TO POINTS

        # TODO: Fields related to edges, 2D cells etc.
        # Notes:
        # u = VN.vtk_to_numpy(self.readerData.GetCellData().GetArray(
        # 'velocity'))
        # u = VN.vtk_to_numpy(self.readerData.GetFieldData().GetArray(
        # 'velocity'))

        # Get all point related data including fields
        pointData = self.readerData.GetPointData()
        # Get number of fields
        num_pointFields =  pointData.GetNumberOfArrays()
        print("Number of VTK fields: \n", num_pointFields)

        # cellData = self.readerData.GetCellData()
        # num_cellFields =  cellData.GetNumberOfArrays()

        # fieldData = self.readerData.GetFieldData()
        # num_fieldData =  fieldData.GetNumberOfArrays()

        # Note: to read all scalars and vectors the next functions can be used:
        # self.readerData.ReadAllScalarsOn()
        # self.readerData.ReadAllVectorsOn()

        # Loop through all fields
        for i in range(num_pointFields):
            # Get array object
            array = pointData.GetArray(i)
            # Get field name
            # fieldName = pointData.GetArrayName(i) # second way of obtaining name
            fieldName = array.GetName()
            print("Reading field: \n", fieldName)
            # size = array.GetSize()
            # Get values of a field
            values = VN.vtk_to_numpy(pointData.GetArray(i))
            # Add new dictionary entry.
            # Key = fieldName,
            # Associated data = values
            self.readerData.field_dict[fieldName] = values

        # In case there is no values in pointData
        # TODO: needs reimplementation
        num_cellFields = 0
        if num_pointFields == 0:
            cellData = self.readerData.GetCellData()
            num_cellFields = cellData.GetNumberOfArrays()

            for i in range(num_cellFields):
                # Get array object
                array = cellData.GetArray(i)
                # Get field name
                fieldName = array.GetName()
                print("Reading field: \n", fieldName)
                # size = array.GetSize()
                # Get values of a field
                values = VN.vtk_to_numpy(cellData.GetArray(i))
                # Add new dictionary entry.
                # Key = fieldName,
                # Associated data = values
                self.readerData.field_dict[fieldName] = values



    def getDataSetType(self):
        """Get VTK file (ASCII) dataset type (e.g. POLYDATA, UNSTRUCTURED_GRID etc.)

        arguments:
            filePath (str) : Full path (including file name) to ASCII TK file.
        """

        # Open file
        vtkFile = open(self.filePath, "r")
        # Set default variable
        dataset_type = "None"
        # Interpolate through lines
        for line in vtkFile:
            if "DATASET" in line:
                l = line.split(' ')
                dataset_type = l[1]
                # Break when dataset has been found
                break

        return dataset_type

class vtkWriteUtilities:
    def __init__(self, fileName="test"):
        self.fileName = fileName

    def writeGridGeometry(self, grid_dict):
        """Write unstructured grid data to .vtk file.

        Arguments:
            grid_dict  (dictionary) : Dictionary holding grid geometry
                                      definition, topology elements definition.
                                      For full description see the output
                                      variable of class "vtkReadUtilities" -->
                                      function "readVtkUnstructuredGrid"

        """
        from pyvtk import VtkData, UnstructuredGrid

        print("vtk_write STARTED")

        points_geo = grid_dict['points_geo']
        obj_0D_list = grid_dict['obj_0D_list']
        obj_1D_list = grid_dict['obj_1D_list']
        obj_2D_list = grid_dict['obj_2D_list']
        obj_3D_list = grid_dict['obj_3D_list']

        print("Number of points: ", len(points_geo))
        print("Number of nodes: ", len(obj_0D_list))
        # Converting 2D points to 3D (addint third coordinate as 0.0)
        if len(points_geo[0] == 2):
            z = np.zeros((len(points_geo),3))
            z[:,:-1] = points_geo
            points_geo = z
        print("Number of edges: ", len(obj_1D_list))
        print("Number of 2D cells: ", len(obj_2D_list))

        num_0D_obj = len(obj_0D_list)
        num_1D_obj = len(obj_1D_list)
        num_2D_obj = len(obj_2D_list)
        num_3D_obj = len(obj_3D_list)

        # Create an execution string, regarding on the types of the elements
        # found within the IDS database
        # Execution string example:
        # structure = UnstructuredGrid(points = points_geo,
        #                              vertex = obj_0D_list,
        #                              line = obj_1D_list,
        #                              polygon = obj_2D_list,
        #                              tetra = obj_3D_list)
        # help: https://github.com/pearu/pyvtk/blob/master/pyvtk/UnstructuredGrid.py
        execution_string = 'self.structure = UnstructuredGrid(points = points_geo'
        # - 0D objects/elements : vertex
        if num_0D_obj > 0:
            execution_string = execution_string + ', vertex = obj_0D_list'
        # - 1D objects/elements : line/edge
        if num_1D_obj > 0:
            execution_string = execution_string + ', line = obj_1D_list'
        # - 2D objects/elements : triangle or quad or polygon
        if num_2D_obj > 0:
            num_points_2D_element = len(obj_2D_list[0])
            if num_points_2D_element == 3:
                execution_string = execution_string + ', triangle = obj_2D_list'
            elif num_points_2D_element == 4:
                execution_string = execution_string + ', quad = obj_2D_list'
            else:
                execution_string = execution_string + ', polygon = obj_2D_list'
        # - 3D objects/elements : tetra or hexahedron
        if num_3D_obj > 0:
            num_points_3D_element = len(obj_3D_list[0])
            if num_points_3D_element == 4:
                execution_string = execution_string + ', tetra = obj_3D_list'
            elif num_points_3D_element == 8:
                execution_string = execution_string + ', hexahedron = obj_3D_list'
        execution_string = execution_string + ')'

        print('Execution string: ', execution_string)
        # - Run the execution string

        # exec('self.structure = UnstructuredGrid(points = points_geo, vertex = obj_0D_list, line = obj_1D_list)')
        exec(execution_string)

        vtk = VtkData(self.structure)
        vtk.tofile(self.fileName+'.vtk','ascii')
        print("vtk_write FINISHED")

    def writeUnstructuredGrid(self, grid_dict, field_dict=None):
        """Write unstructured grid geometry and field data to .vtk file.

        Arguments:
            grid_dict  (dictionary) : Dictionary holding grid geometry
                                      definition, topology elements definition.
                                      For full description see the output
                                      variable of class "vtkReadUtilities" -->
                                      function "readVtkUnstructuredGrid".
            field_dict (dictionary) : Dictionary holding all field data.
                                      For full description see the output
                                      variable of class "vtkReadUtilities" -->
                                      function "readVtkDataFields".


        """

        # Set vtk objects for unstructured grid, points and cell arrays
        ug = vtk.vtkUnstructuredGrid()
        points = vtk.vtkPoints()
        cells = vtk.vtkCellArray()

        # Set points
        for point in grid_dict["points_geo"]:

            if len(point) == 2:
                # For 2D set third coordinate as 0.0
                point = np.append(point, 0.0)
                points.InsertNextPoint(point)
            else:
                points.InsertNextPoint(point)

        print("Setting data fields")

        if field_dict != None:

            # Set points field data
            for fkey in field_dict.keys():
                # Set array of float/double values
                scalars_points = vtk.vtkFloatArray()
                # scalars_points = vtk.vtkDoubleArray()
                scalars_points.SetNumberOfComponents(1)
                # Set field name
                scalars_points.SetName(fkey)
                # Get number of scalars
                num_scalars = len(field_dict[fkey])

                print("- " + fkey + ": " + str(num_scalars) + " scalars")
                # Set scalars
                for i in range(num_scalars):
                    scalars_points.InsertNextTuple1(field_dict[fkey][i])
                # Set new array to unstructured grid
                ug.GetPointData().AddArray(scalars_points)

                print("  Range: ", scalars_points.GetRange())


        # TODO: for other type of cells (lines etc.)

        # cell_type =

        i = 0
        for cell in grid_dict["obj_2D_list"]:
            cells.InsertNextCell(4)
            cells.InsertCellPoint(grid_dict["obj_2D_list"][i][0])
            cells.InsertCellPoint(grid_dict["obj_2D_list"][i][1])
            cells.InsertCellPoint(grid_dict["obj_2D_list"][i][2])
            cells.InsertCellPoint(grid_dict["obj_2D_list"][i][3])

            i += 1

        ug.SetPoints(points)
        ug.SetCells(9, cells)

        writer = vtk.vtkDataSetWriter()
        writer.SetFileName(self.fileName+'.vtk')
        writer.SetInputData(ug)
        writer.SetHeader("Header text")
        writer.Write()

        # vtkDataSetWriter writes values using comma as decimal separator.
        # This is then unreadable by ParaView!!!!
        # There was no option found to change that! Open the file and change
        # all ',' to '.'
        self.vtkFileComma2Period()

    def vtkFileComma2Period(self):
        """In file change all commas (",") to periods (".").
        """

        import fileinput

        with fileinput.FileInput(self.fileName+'.vtk',
                                 inplace=True,
                                 backup='.bak') as file:
            for line in file:
                print(line.replace(',', '.'), end='')
