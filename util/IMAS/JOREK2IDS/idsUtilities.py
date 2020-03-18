#!/usr/bin/env python

#  Name : jorekHDF5toIDS.py
#
#  Description :
#           A collection of routines for handling IDSs.
#
#           Requirements:
#               - IMAS (python library!)
#
#  Author :
#         Dejan Penko
#  E-mail :
#         dejan.penko@lecad.fs.uni-lj.si
#
#****************************************************
#     Copyright(c) 2020- D. Penko

ENABLED = True
try:
    import imas
except ImportError as e:
    # No IMAS module installed
    ENABLED = False

class basicIDS(object):
    """This class provides basic routines for opening and closing the IMAS
    databases - IDSs.
    """

    def __init__(self, shot, run, user, device, version='3'):
        """
        Arguments:
            shot    (int): Shot number of the IDS case.
            run     (int): Run number of the IDS case.
            user    (str): IMAS database owner user name.
            device  (str): Device name.
            version (str): IMAS major version.
            grid_ggdName (str): Name of grid_ggd from which data was retrieved.
        """
        self.shot = shot
        self.run = run
        self.user = user
        self.device = device
        self.version = version
        self.state = False
        self.ggdName = ''
        # IDS data entry for calling put
        self.ids = None
        # GGD for allocating and storing mesh data
        self.grid_ggd = None
        self.ggd = None
        self.imas_obj = None

    def openIMASdatabase(self):
        """Open IMAS database.
        """

        print('Open IMAS database START')
        self.state = True

        if not ENABLED:
            self.state = False
            return

        try:
            self.imas_obj = imas.ids(self.shot, self.run, self.shot, self.run)
            self.imas_obj.open_env(self.user, self.device, self.version)
            if self.imas_obj.isConnected():
                print('Creation of data entry OK!')
            else:
                print('Creation of data entry FAILED!')
                self.state = False
        except Exception as e:
            print('Failed to open IMAS database. The specified IMAS ' + \
                  'database does not exist!')
            self.state = False

        print('Open IMAS database FINISHED')

    def createIMASdatabase(self):
        """Created IMAS database.
        """
        print('Create IMAS database START')
        self.state = True

        if not ENABLED:
            self.state = False
            return
        try:
            self.imas_obj = imas.ids(self.shot, self.run, self.shot, self.run)
            self.imas_obj.create_env(self.user, self.device, self.version)
            if self.imas_obj.isConnected():
                print('Creation of data entry Ok!')
                self.state = True
            else:
                print('Creation of data entry FAILED!')
                self.state = False
        except Exception as e:
            print('Failed to create IMAS database. Possible problems with IMAS'
                  ' module!')
            self.state = False
        print('Create IMAS database FINISHED')

    def idsClose(self):
        """Close IDS database"""
        self.imas_obj.close()


class readIDS(basicIDS):
    """This class provides basic routines for opening and closing the IMAS
    databases - IDSs.
    """

    def __init__(self, shot, run, user, device, version='3'):
        super(readIDS, self).__init__(shot, run, user, device, version)
        if ENABLED:
            self.openIMASdatabase()

    def getGGD(self, IDSName):
        """Get GGD from the "IDSName" IDS (e.g. "edge_profiles, wall etc.).
        Returns the first (0) slice of ggd. For non-dynamic there is only
        one slice. So by default the first
        slice is returned.

        Arguments:
            IDSName  (str)      : Name of the IDS (e.g. edge_profiles)

        Returns:
            grid_ggd (ids.grid_ggd) : GRID GGD object.

        """

        self.ggdName = IDSName

        print('Get IDS GGD START')

        # Define self.grid_ggd
        if IDSName == 'edge-profiles':
            # Get the data and set 'Shortcut' variable to imas_obj.edge_profiles
            # data tree node
            imas_obj.edge_profiles.get()
            ids = imas_obj.edge_profiles
            # 'Shortcut' variable to ...grid_ggd[0] node
            self.grid_ggd = ids.grid_ggd
        elif IDSName == 'wall':
            # Get the data and set 'Shortcut' variable to imas_obj.wall
            # data tree node
            imas_obj.wall.get()
            ids = imas_obj.wall

            # 'Shortcut' variable to ...grid_ggd node
            self.grid_ggd = ids.description_ggd[0].grid_ggd
            self.ggd = ids.description_ggd[0].ggd
        else:
            print('The specified IDS either is not supported or it does ' + \
                  'not exist')
            self.grid_ggd = None

        # Return False if self.grid_ggd was not set
        if self.grid_ggd == None:
            return False

        print('Get IDS GGD FINISHED')
        return True

    def getIDSMesh(self, n=0):
        """Read IDS, get and set the necessary grid data (points geometry,
        connectivity array for 1D and 2D objects).

        Arguments:
            n (int): The index of the slice in GGD
        """
        print("getMesh START")
        grid_ggd = self.grid_ggd[n]

        points_geo = []
        list0D = []    # vertices
        list1D = []    # edges
        list2D = []    # 2D cells /faces
        list3D = []    # 3D volumes

        # Get name of the mesh
        meshName = str(ggd.identifier.name)

        # - Get number of 0D objects - grid nodes
        # Note: This list contains points for all objects. First nodes, then
        # edges, then 2D cells/faces and lastly 3D volumes.
        # Even if the objects share the same point, for both
        # objects is the point defined separately.
        # Get the max dimension of the grid.
        DIM = len(ggd.space[0].objects_per_dimension)
        # - Get number of 0D objects - nodes
        if DIM >= 1:
            # - 'Shortcut' variable to ...objects_per_dimension[0] node
            ids_dim_0D = ggd.space[0].objects_per_dimension[0]
            num_obj_0D_all = len(ids_dim_0D.object)
        else:
            num_obj_0D_all = 0
        if DIM >= 2:
            # - 'Shortcut' variable to ...objects_per_dimension[1] node
            ids_dim_1D = ggd.space[0].objects_per_dimension[1]
            # - Get number of 1D objects - edges
            num_obj_1D_all = len(ids_dim_1D.object)
        else:
            num_obj_1D_all = 0
        if DIM >= 3:
            # - 'Shortcut' variable to ...objects_per_dimension[2] node
            ids_dim_2D = ggd.space[0].objects_per_dimension[2]
            # - Get number of 2D objects - cells/faces
            num_obj_2D_all = len(ids_dim_2D.object)
        else:
            num_obj_2D_all = 0
        if DIM >= 4:
            # - 'Shortcut' variable to ...objects_per_dimension[3] node
            ids_dim_3D = ggd.space[0].objects_per_dimension[3]
            # - Get number of 3D objects - volumes
            num_obj_3D_all = len(ids_dim_3D.object)
        else:
            num_obj_3D_all = 0

        print("num_obj_0D_all: ", num_obj_0D_all)
        print("num_obj_1D_all: ", num_obj_1D_all)
        print("num_obj_2D_all: ", num_obj_2D_all)
        print("num_obj_3D_all: ", num_obj_3D_all)

        # Get geometry (coordinates) of each point and
        # get connectivity array (node list) of the 0D objects
        for i in range(num_obj_0D_all):
            points_geo.append(ids_dim_0D.object[i].geometry)
            list0D.append(int(ids_dim_0D.object[i].nodes[0]))
            # Note: Transforming index to python notation (vtk uses that)
            #       Python_index = Fortran_index - 1
        # Get connectivity array (node list) of the 1D objects
        if num_obj_1D_all > 0 and ids_dim_1D.object[0].nodes[0] != 0:
            object = ids_dim_1D.object
            for i in range(num_obj_1D_all):
                # Create a list of all objects the 2D object is comprised of
                el1D = [int(object[i].nodes[j]) for j in range(2)]
                # Add the object array of indices to the connectivity array
                list1D.append(el1D)
        # Get connectivity array (node list) of the 2D objects
        # - Get the shape of the 2D objects (triangle, quad...). it is presumed
        #   that all 2D objects are of the same shape.
        if num_obj_2D_all > 0 and ids_dim_2D.object[0].nodes[0] != 0:
            object = ids_dim_2D.object
            m = len(object[0].nodes)
            for i in range(num_obj_2D_all):
                # Create a list of all objects the 2D object is comprised of
                el2D = [int(object[i].nodes[j]) for j in range(m)]
                # Add the object array of indices to the connectivity array
                list2D.append(el2D)

        # Get connectivity array (node list) of the 3D objects
        # - Get the shape of the 3D objects (tetra, hexahedron...). it is presumed
        #   that all 3D objects are of the same shape.
        if num_obj_3D_all > 0 and ids_dim_3D.object[0].nodes[0] != 0:
            object = ids_dim_3D.object
            m = len(object[0].nodes)
            for i in range(num_obj_3D_all):
                # Create a list of all objects the 3D object is comprised of
                el3D = [int(object[i].nodes[j]) for j in range(m)]
                # Add the object array of indices to the connectivity array
                list3D.append(el3D)

        # print("First list0D: ", list0D[0])
        # print("First obj_1D: ", list1D[0])
        # print("First obj_2D: ", list2D[0])

        print("getMesh from slice %d FINISHED" % n)

        return points_geo, list0D, list1D, list2D, list3D, meshName

    def createMesh(self):
        if self.grid_ggd is None:
            print('No GGD has been retrieved!')
            return

        # Get number of slices
        N = len(self.grid_ggd)
        meshes = []

        for i in range(N):
            nodes, El0D, El1D, El2D, El3D, meshName = self.getIDSMesh(i)

            # Create new mesh
            if meshName == '':
                meshName = 'IDS-%d-%d-%s-slice-%d' % (self.shot, self.run,
                                                      self.ggdName, i)
            meshes.append(meshName)

            SMESH_utils.createSMesh(nodes, El0D, El1D, El2D, El3D, meshName)
        return True, meshes

    @staticmethod
    def checkEqdskIDS(ids_equilibrium, ids_wall):
        """Check if IDS has valid equilibrium (eqdsk) data.
        Note: Checks only the critical ones.
        """

        if len(ids_equilibrium.time_slice[0].boundary_separatrix.outline.r) != 1:
            err_msg = 'Faulty RDIM found in IDS.'
            return False, err_msg
        if len(ids_equilibrium.time_slice[0].boundary_separatrix.outline.z) != 1:
            err_msg = 'Faulty ZDIM found in IDS.'
            return False, err_msg
        if len(ids_equilibrium.vacuum_toroidal_field.b0) != 1:
            err_msg = 'Faulty BCENTR found in IDS.'
            return False, err_msg

        return True, ''

    def getEqdsk(self):
        """Read Equilibrium and Wall IDS to get G-EQDSK file format relevant
        data.
        """
        import eqdsk
        ids_equilibrium = self.imas_obj.equilibrium
        ids_equilibrium.get()
        ids_wall = self.imas_obj.wall
        ids_wall.get()

        ok, err_msg = self.checkEqdskIDS(ids_equilibrium, ids_wall)
        if not ok:
            self.displayMessage('Error', err_msg)
            return False

        self.x = eqdsk.EQDSK()

        # Get eqdsk data from Equilibrium IDS
        self.x.setHEADER(ids_equilibrium.ids_properties.comment )
        self.x.setRDIM(ids_equilibrium.time_slice[0].boundary_separatrix.outline.r)
        self.x.setZDIM(ids_equilibrium.time_slice[0].boundary_separatrix.outline.z)
        self.x.setRCENTR(ids_equilibrium.vacuum_toroidal_field.r0)
        self.x.setRLEFT(ids_equilibrium.time_slice[0].boundary.minor_radius)
        self.x.setZMID(ids_equilibrium.time_slice[0].boundary.geometric_axis.z)
        self.x.setRMAXIS(ids_equilibrium.time_slice[0].global_quantities.magnetic_axis.r)
        self.x.setZMAXIS(ids_equilibrium.time_slice[0].global_quantities.magnetic_axis.z)
        self.x.setSIMAG(ids_equilibrium.time_slice[0].global_quantities.psi_axis)
        self.x.setSIBRY(ids_equilibrium.time_slice[0].global_quantities.psi_boundary)
        self.x.setBCENTR(ids_equilibrium.vacuum_toroidal_field.b0)
        self.x.setCURRENT(ids_equilibrium.time_slice[0].global_quantities.ip )
        self.x.setFPOL(ids_equilibrium.time_slice[0].profiles_1d.f)
        self.x.setPRES(ids_equilibrium.time_slice[0].profiles_1d.pressure )
        self.x.setFFPRIM(ids_equilibrium.time_slice[0].profiles_1d.f_df_dpsi)
        self.x.setPPRIME(ids_equilibrium.time_slice[0].profiles_1d.dpressure_dpsi)
        self.x.setPSIRZ(ids_equilibrium.time_slice[0].profiles_2d[0].psi)
        self.x.setQPSI(ids_equilibrium.time_slice[0].profiles_1d.q)
        self.x.setRBBBS(ids_equilibrium.time_slice[0].boundary.outline.r)
        self.x.setZBBBS(ids_equilibrium.time_slice[0].boundary.outline.z)
        self.x.setNBBBS(len(self.x.getRBBBS()))
        self.x.setNW(len(self.x.getPRES()))
        self.x.setNH(len(self.x.getPSIRZ()))

        # Get eqdsk related data from Wall IDS
        self.x.setRLIM(ids_wall.description_2d[0].limiter.unit[0].outline.r)
        self.x.setZLIM(ids_wall.description_2d[0].limiter.unit[0].outline.z)
        self.x.setLIMITR(len(self.x.getRLIM()))

        self.x.successfullRead = True

        return True

    def createEqdsk(self, path):
        """Create G-EQDSK format file from data obtained from Equilibrium and
        Wall IDSs.
        """
        # import eqdsk
        # x = eqdsk.EQDSK()
        eqdskFile = str(path)
        eqdskName = eqdskFile.rsplit('/', 1)[-1] # Get the name of the file only
        self.x.setName(eqdskName)
        # The values can be changed and you can generate a EQDSK G file with the
        # following function
        # text = self.x.generateText()
        # This generates the string so you have to manually save the text in a
        # file
        with open(eqdskFile, 'w') as f:
            f.write(self.x.generateText())
        return True

class writeIDS(basicIDS):
    def __init__(self, shot, run, user, device, version='3'):
        super(writeIDS, self).__init__(shot, run, user, device, version)
        if ENABLED:
            self.createIMASdatabase()

    def createGridGGD(self, IDSName, numSlices):
        """Open an IDS.
        """
        self.grid_ggdName = IDSName
        print('Create IDS GGD START')

        if IDSName == "wall":
            self.ids, self.grid_ggd = self.createGridGGDWall(self.imas_obj,
                                                    numSlices)
        elif IDSName == 'edge_profiles':
            self.ids, self.grid_ggd = self.createGridGGDEdgeProfiles(self.imas_obj,
                                                            numSlices)

        elif IDSName == 'mhd':
            self.ids, self.grid_ggd = self.createGridGGDMhd(self.imas_obj, numSlices)

        else:
            print('The specified IDS either is not supported or it does not '
                  'exist')
            return False
        print('Create IDS GGD FINISHED')
        return True

    def setBaseEquilibrium(self):
        """Set fundamental data for Equilibrium IDS.
        """
        time = 0.0
        ids_equilibrium = self.imas_obj.equilibrium
        ids_equilibrium.putNonTimed()
        ids_equilibrium.ids_properties.homogeneous_time = 1
        ids_equilibrium.time = [time]
        return ids_equilibrium

    def setBaseWall(self):
        """Set fundamental data for Wall IDS.
        """
        time = 0.0
        ids_wall = self.imas_obj.wall
        ids_wall.putNonTimed()
        ids_wall.ids_properties.homogeneous_time = 1
        ids_wall.time.resize(1)
        ids_wall.time = [time]
        return ids_wall

    def createDatasetDescriptionIDS(self, user, source='', comment=''):
        """ Set data to be written to the dataset_description IDS,
        designed to hold this IDS basic information (IMAS version used
        for writing the database and similar).
        """
        import time
        import os
        # 'Shortcut' variable to imas_obj.dataset_description IDS
        desc = self.imas_obj.dataset_description
        # Set IMAS database description and base information
        props = desc.ids_properties
        # props.comment = ""
        props.homogeneous_time = 1
        props.source = source
        props.comment = comment
        props.provider = user
        props.creation_date = time.strftime("%d-%m-%Y")
        desc.imas_version = os.environ['IMAS_VERSION']
        desc.dd_version = os.environ['UAL_VERSION']
        desc.time.resize(1)
        desc.time = [time]

        # Write the database description information
        desc.putSlice()

        return True

    @staticmethod
    def createGridGGDEdgeProfiles(imas_obj, numSlices):
        time = 0.0
        ids = imas_obj.edge_profiles
        ids.putNonTimed()
        ids.ids_properties.homogeneous_time = 1
        ids.time.resize(1)
        # ids.time = [time]
        ids.time = [0]
        # Resize the grid_ggd[:] node
        ids.grid_ggd.resize(numSlices)

        return ids, ids.grid_ggd

    @staticmethod
    def createGridGGDWall(imas_obj, numSlices):
        time = 0.0
        ids = imas_obj.wall
        ids.putNonTimed()
        ids.ids_properties.homogeneous_time = 1
        ids.time.resize(1)
        ids.time = [time]
        # Resize the grid_ggd[:] node
        ids.description_ggd.resize(1)
        ids.description_ggd[0].grid_ggd.resize(numSlices)
        return ids, ids.description_ggd[0].grid_ggd

    @staticmethod
    def createGridGGDMhd(imas_obj, numSlices):
        time = 0.0
        ids = imas_obj.mhd
        # ids.putNonTimed()
        ids.ids_properties.homogeneous_time = 1
        ids.time.resize(1)
        ids.time = [time]
        # Resize tree node (allocate)
        ids.grid_ggd.resize(numSlices)
        return ids, ids.grid_ggd

    def writeMeshToSlice(self, points_geo, obj_0D_list, obj_1D_list,
                         obj_2D_list, obj_3D_list, n = 0, label=''):
        """Prepare data to be written to GGD, in the n-th slice.
        """

        if self.grid_ggd[n] is None:
            print('No IDS database created! Aborting')
            return False

        grid_ggd = self.grid_ggd[n]

        # Set size of the grid_ggd[:].space[:] node
        grid_ggd.space.resize(1)
        grid_ggd.space[0].objects_per_dimension.resize(4)

        # Set the coordinate system information
        # - Set number of coordinate types
        num_coordtype = 3
        # - Set coordinates_type
        grid_ggd.space[0].coordinates_type.resize(num_coordtype)
        # - Fill coordinates_type
        grid_ggd.space[0].coordinates_type = [1, 2, 3]

        # Additional information on the grid
        # TODO

        grid_description = "Mesh %s" % label
        grid_ggd.identifier.description = grid_description
        grid_ggd.identifier.name = label
        grid_ggd.identifier.index = 1

        # Writing NODES
        num_points = len(points_geo)

        ids_space = grid_ggd.space[0]
        ids_space.objects_per_dimension[0].object.resize(num_points)
        ids_dim_0D = ids_space.objects_per_dimension[0]
        ids_dim_1D = ids_space.objects_per_dimension[1]
        ids_dim_2D = ids_space.objects_per_dimension[2]
        ids_dim_3D = ids_space.objects_per_dimension[3]

        # # Writing 0D elements
        for i in range(num_points): # it should be num_points == len(num_obj_0D_all)
            ids_dim_0D.object[i].geometry.resize(num_coordtype)
            ids_dim_0D.object[i].nodes.resize(1)
            ids_dim_0D.object[i].geometry = points_geo[i]
            ids_dim_0D.object[i].nodes[0] = i + 1

        num_obj_0D_all = len(obj_0D_list)
        num_obj_1D_all = len(obj_1D_list)
        num_obj_2D_all = len(obj_2D_list)
        num_obj_3D_all = len(obj_3D_list)

        # ids_dim_0D.object.resize(num_obj_0D_all)

        # Writing 1D objects
        if num_obj_1D_all > 0:
            ids_dim_1D.object.resize(num_obj_1D_all)
            for i in range(num_obj_1D_all):
                object = ids_dim_1D.object[i]
                object.nodes.resize(2)
                for j in range(2):
                    if obj_1D_list[i][j] == 0:
                        print("WARNING! Index 0 found while writing 1D objects!")
                    object.nodes[j] = obj_1D_list[i][j]
        else:
            # Write a dummy object, holding value 0 in nodes array.
            # This is required in case the mesh does not contain 1D objects
            # but does contain 2D or 3D object. The absence of this dummy
            # would result in objects_per_dimension[1] not being skipped but
            # as an 'end to write to objects_per_dimension structure' indicator!
            ids_dim_1D.object.resize(1)
            object = ids_dim_1D.object[0]
            object.nodes.resize(1)
            object.nodes[0] = 0

        # Writing 2D objects
        if num_obj_2D_all > 0:
            ids_dim_2D.object.resize(num_obj_2D_all)
            m = len(obj_2D_list[0])
            for i in range(num_obj_2D_all):
                object = ids_dim_2D.object[i]
                object.nodes.resize(m)
                for j in range(m):
                    if obj_2D_list[i][j] == 0:
                        print("WARNING! Index 0 found while writing 2D objects!")
                    # Set each objects nodes data
                    object.nodes[j] = obj_2D_list[i][j]
        else:
            # Write a dummy object, holding value 0 in nodes array.
            # This is required in case the mesh does not contain 1D objects
            # but does contain 3D object. The absence of this dummy
            # would result in objects_per_dimension[2] not being skipped but
            # as an 'end to write to objects_per_dimension structure' indicator!
            ids_dim_2D.object.resize(1)
            object = ids_dim_2D.object[0]
            object.nodes.resize(1)
            object.nodes[0] = 0

        # Writing 3D objects
        if num_obj_3D_all > 0:
            ids_dim_3D.object.resize(num_obj_3D_all)
            m = len(obj_3D_list[0])
            for i in range(num_obj_3D_all):
                object = ids_dim_3D.object[i]
                object.nodes.resize(m)
                for j in range(m):
                    # Set each objects nodes data
                    object.nodes[j] = obj_3D_list[i][j]


        # WRITING GRID SUBSETS
        # Setting initial number of grid subsets = 1 (for all points in the domain)
        num_grid_subset = 1
        if num_obj_0D_all > 0:
            num_grid_subset += 1
        if num_obj_1D_all > 0:
            num_grid_subset += 1
        if num_obj_2D_all > 0:
            num_grid_subset += 1
        # if num_obj_3D_all > 0:
        #     num_grid_subset += 1

        print("Number of grid_subsets to be made: ", str(num_grid_subset))
        grid_ggd.grid_subset.resize(num_grid_subset)
        gs_index = 0

        # 0D objects/points
        # Note: num_points refer to points, num_obj_0D_all refers to vertices
        # (VTK separates points and vertices)
        if num_points > 0:
            gs = grid_ggd.grid_subset[gs_index]
            gs.identifier.name = "Nodes "
            gs.identifier.index = gs_index + 1 # Fortran notation
            gs.identifier.description = """All points/nodes/vertices/0D
                                      objects in the domain."""
            gs.dimension = 1
            gs.element.resize(num_points)
            for j in range(num_points):
                gs.element[j].object.resize(1)
                # Write in Fortran notation (!)
                gs.element[j].object[0].space = 1
                gs.element[j].object[0].index = j + 1 # Fortran notation
                gs.element[j].object[0].dimension = 1

            gs_index += 1

        # 1D objects
        if num_obj_1D_all > 0:
            gs = grid_ggd.grid_subset[gs_index]
            gs.identifier.name = "Edges"
            gs.identifier.index = gs_index + 1 # Fortran notation
            gs.identifier.description = "All edges/1D objects in the domain."
            gs.dimension = 2
            gs.element.resize(num_obj_1D_all)
            for j in range(num_obj_1D_all):
                gs.element[j].object.resize(1)
                # Write in Fortran notation (!)
                gs.element[j].object[0].space = 1
                gs.element[j].object[0].index = j + 1 # Fortran notation
                gs.element[j].object[0].dimension = 2

            gs_index += 1

        # 2D objects
        if num_obj_2D_all > 0:
            gs = grid_ggd.grid_subset[gs_index]
            gs.identifier.name = "2D cells"
            gs.identifier.index = gs_index + 1 # Fortran notation
            gs.identifier.description = "All 2D cells/2D objects in the domain."
            gs.dimension = 3
            gs.element.resize(num_obj_2D_all)
            for j in range(num_obj_2D_all):
                gs.element[j].object.resize(1)
                # Write in Fortran notation (!)
                gs.element[j].object[0].space = 1
                gs.element[j].object[0].index = j + 1 # Fortran notation
                gs.element[j].object[0].dimension = 3

            gs_index += 1

        return True

    def ggdWriteQuantityArray(self, IDSQuantityPath, array, gsInd, slice=0):
        """Write 1D scalar array, containing quantity values which correspond to
        a grid_subset.

        Arguments:
            IDSQuantityPath  (obj)      : Path to the IDS quantity tree node e.g.
                                  ...ggd[0].electrons[0].temperature[0]
            array (array)    : 1D scalar array which contains quantity
                                  values, corresponding to grid subset with
                                  index gs_ind
            gsInd   (int)      : Index of grid subset, to which the 1D scalar
                                  array corresponds to.
            slice:  (int)      : GGD slice.
        """


        IDSQuantityPath.grid_index = 1
        IDSQuantityPath.grid_subset_index = gsInd
        IDSQuantityPath.values = array

        return True

    def writeMesh(self, meshList):
        """Writes the mesh to IDS
        """
        for i in range(len(meshList)):
            if meshList[i] == '':
                print('No mesh name provided')
                continue
            # Write mesh to GGD, slice i
            ok = self.writeMeshToSlice(meshList[i], i)
            if not ok:
                print('Failed to write mesh %s to slice %d' %
                      meshList[i], i)

        # Try putting data to ids
        try:
            self.ids.put()
        except Exception as e:
            print('Error when trying to write data to IDS!')
            self.idsClose()
            return False
        print('Successfully written data to IDS!')
        self.idsClose()
        return True

    def writeEqdsk(self, eqdsk):
        """Write the data obtained from the G-EQDSK format file to Equilibrium
         and Wall IDS.

        Arguments:
            eqdsk (obj) : Eqdsk.py script object.
        """
        import time

        # Set Equilibrium IDS fundamental data
        ids_equilibrium = self.imas_obj.equilibrium # Equilibrium IDS object

        # Set Equilibrium IDS properties
        ids_equilibrium.ids_properties.provider = self.user
        ids_equilibrium.ids_properties.source = 'EQDSK2IDS SMITER utility. ' \
            'This Equilibrium IDS holds the data of the equilibrium (G-EQDSK) '\
            ' file: ' + eqdsk.getName() + '.'
        ids_equilibrium.ids_properties.comment = eqdsk.getHEADER()

        ids_equilibrium.ids_properties.creation_date = time.strftime("%d-%m-%Y")

        # Resize time slice structure
        ids_equilibrium.time_slice.resize(1)
        # Set eqdsk data to Equilibrium IDS
        ids_equilibrium.time_slice[0].boundary_separatrix.outline.r = [eqdsk.getRDIM()]
        ids_equilibrium.time_slice[0].boundary_separatrix.outline.z = [eqdsk.getZDIM()]
        ids_equilibrium.vacuum_toroidal_field.r0 = eqdsk.getRCENTR()
        ids_equilibrium.time_slice[0].boundary.minor_radius = eqdsk.getRLEFT()
        ids_equilibrium.time_slice[0].boundary.geometric_axis.z = eqdsk.getZMID()
        ids_equilibrium.time_slice[0].global_quantities.magnetic_axis.r = eqdsk.getRMAXIS()
        ids_equilibrium.time_slice[0].global_quantities.magnetic_axis.z = eqdsk.getZMAXIS()
        ids_equilibrium.time_slice[0].global_quantities.psi_axis = eqdsk.getSIMAG()
        ids_equilibrium.time_slice[0].global_quantities.psi_boundary = eqdsk.getSIBRY()
        ids_equilibrium.vacuum_toroidal_field.b0 = [eqdsk.getBCENTR()]
        ids_equilibrium.time_slice[0].global_quantities.ip  = eqdsk.getCURRENT()
        ids_equilibrium.time_slice[0].profiles_1d.f = eqdsk.getFPOL()
        ids_equilibrium.time_slice[0].profiles_1d.pressure  = eqdsk.getPRES()
        ids_equilibrium.time_slice[0].profiles_1d.f_df_dpsi = eqdsk.getFFPRIM()
        ids_equilibrium.time_slice[0].profiles_1d.dpressure_dpsi = eqdsk.getPPRIME()
        ids_equilibrium.time_slice[0].profiles_2d.resize(1)
        ids_equilibrium.time_slice[0].profiles_2d[0].psi = eqdsk.getPSIRZ()
        ids_equilibrium.time_slice[0].profiles_1d.q = eqdsk.getQPSI()
        # NBBBS = len(RBBBS)
        # NW = len(PRES)
        # NH = len(PSIRZ)
        ids_equilibrium.time_slice[0].boundary.outline.r = eqdsk.getRBBBS()
        ids_equilibrium.time_slice[0].boundary.outline.z = eqdsk.getZBBBS()

        # Try putting data to Equilibrium IDS
        try:
            ids_equilibrium.put()
        except Exception as e:
            print('Error when trying to write data to Equilibrium IDS!')
            self.idsClose()
            return False

        # Set Wall IDS fundamental data
        ids_wall = self.imas_obj.wall # Wall IDS object
        # Set eqdsk related data to Wall IDS
        ids_wall.description_2d.resize(1)
        ids_wall.description_2d[0].limiter.unit.resize(1)
        ids_wall.description_2d[0].limiter.unit[0].outline.r = eqdsk.getRLIM()
        ids_wall.description_2d[0].limiter.unit[0].outline.z = eqdsk.getZLIM()
        # LIMITR = len(RLIM)

        # Try putting data to Wall IDS
        try:
            ids_wall.put()
        except Exception as e:
            print('Error when trying to write data to Wall IDS!')
            self.idsClose()
            return False

        print('Successfully written data to IDS!')
        self.idsClose()
        return True
