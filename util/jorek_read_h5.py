#!/usr/bin/env python
"""
jorek_read_h5.py

Read JOREK hdf5 restart files and define functions for interpolation
on the exported fields.

Before reading, set:
    - Variable numbers
Before interpolating, set:
    - Without_n0_mode (default false)

Created by Daan van Vugt on 2017-01-18
"""

from __future__ import print_function
import h5py
import numpy as np

prec=np.float32

"""
Class encapsulating some read data from a restart file

input arguments:
    Variables: a list of all variables you want
    i_plane: which plane to interpolate at. If -1 make a 3D field.
"""
class fields(object):
    var_names = ["psi", "u", "j", "w", "rho", "T", "v_par"] # rest is ambiguous
    def read(self, filename, variables=''):
        if (isinstance(variables, str)):
            self.vars = [int(x)-1 for x in variables.split(',')]
        else:
            raise "Error: expected a string of variables"

        with h5py.File(filename, 'r') as hf:
            self.n_var        = hf.get('n_var')[0]
            self.n_period     = hf.get('n_period')[0]
            self.n_tor        = hf.get('n_tor')[0]
            self.n_vertex_max = hf.get('n_vertex_max')[0]
            self.n_elements   = hf.get('n_elements')[0]
            self.vertex       = np.array(hf.get('vertex'), dtype=np.int32)-1
            self.x            = np.array(hf.get('x'), dtype=prec)
            self.size         = np.array(hf.get('size'), dtype=prec)
            if (len(self.vars) > 0):
                self.values   = np.array(hf.get('values')[self.vars,:,:], dtype=prec)

    """
    Create VTK objects from points and connectivity matrix
    input:
        n_sub: number of subdivisions per element
        n_plane: number of planes in toroidal direction
        phi: range (2 elements) or single value of phi
        without_n0_mode: do not include n0 mode if true

    returns:
        vtkUnstructuredGrid
    """
    def to_vtk(self, n_sub=4, phi=[0,2*np.pi], n_plane=16, without_n0_mode=False):
        import vtk
        from vtk.util import numpy_support as npvtk

        # If we do not have a range in phi make only one plane
        periodic=False
        if (isinstance(phi,int)):
            n_plane = 1
            phis = np.asarray([phi])
        elif (n_plane == 1):
            phis = np.asarray([phi[0]])
        else:
            periodic = (np.mod(phi[0]-phi[1],2*np.pi) < 1e-9)
            phis = np.linspace(phi[0],phi[1],num=n_plane,endpoint=not periodic)

        (xyz, ien) = create_grid(self.x, self.vertex, self.size, self.n_elements,
                                 n_sub, phis, n_plane, periodic)

        pcoords = npvtk.numpy_to_vtk(xyz, deep=True, array_type=vtk.VTK_FLOAT)
        points = vtk.vtkPoints()
        points.SetData(pcoords)

        ug = vtk.vtkUnstructuredGrid()
        ug.SetPoints(points)

        HZ = toroidal_basis(self.n_tor, self.n_period, phis, without_n0_mode)
        for i in self.vars:
            val = npvtk.numpy_to_vtk(interp_scalars_3D(self.values[i:i+1,:,:,:],
                                   self.vertex, self.size, n_sub, HZ).reshape(-1),
                                   deep=True, array_type=vtk.VTK_FLOAT)

            val.SetName(self.var_names[i])
            ug.GetPointData().SetScalars(val)


        if (n_plane > 1):
            etype = vtk.VTK_HEXAHEDRON
        else: # or stay 2D
            etype = vtk.VTK_QUAD

        cells = vtk.vtkCellArray()
        cells.SetCells(ien.shape[0], npvtk.numpy_to_vtk(ien, deep=True, array_type=vtk.VTK_ID_TYPE))

        ug.SetCells(etype, cells)
        return ug
    

"""
Calculate RZ positions of all points
return x[element,is,it,var] (where var is 0->R or 1->Z)
"""
def grid_2D(x, vertex, size, n_sub):
    # Calculate RZ for all of the elements (dimension 0) for each of the s
    # positions (dimension 1) for each of the t positions (dimension 2)
    # Multiply x[var, order, node[vertex, element]] on the last two dimensions
    # with size[order, vertex, element]*bf[order, vertex, s, t]
    # See http://stackoverflow.com/questions/26089893/understanding-numpys-einsum
    #return np.einsum('lijk,ijk,ijmn->kmnl', x[:,:,vertex], size, bf(n_sub))
    # Code below is ~5x faster or so! try again when einsum supports optimize=True

    # First create a temporary array holding: x[order, vertex, element, var]
    tmp = np.zeros((x.shape[0], x.shape[1],vertex.shape[0],vertex.shape[1]), dtype=prec)
    # Fill it with the right x
    for i in range(vertex.shape[0]): # small loop over vertices (hardcode 4 here?)
        tmp[:,:,i,:] = x[:,:,vertex[i,:]]
    # multiply by size[order, vertex, element]
    tmp[0,:,:,:] *= size
    tmp[1,:,:,:] *= size
    # Create output array
    out = np.zeros((vertex.shape[1],n_sub,n_sub,2), dtype=prec)
    out[:,:,:,0] = np.tensordot(tmp[0,:,:,:], bf(n_sub), axes=((0,1),(0,1)))
    out[:,:,:,1] = np.tensordot(tmp[1,:,:,:], bf(n_sub), axes=((0,1),(0,1)))
    return out


"""
Calculate toroidal basis functions
"""
def toroidal_basis(n_tor, n_period, phis, without_n0_mode):
    # Setup toroidal coefficients for each plane and toroidal harmonic
    HZ = np.zeros((n_tor,len(phis)))
    for i in range(n_tor):
        mode = np.floor((i+1)/2)*n_period
        if (i == 0):
            if (not without_n0_mode):
                HZ[i,:] = 1
        elif (i % 2 == 0):
            HZ[i,:] = np.sin(mode*phis)
        elif (i % 2 == 1):
            HZ[i,:] = np.cos(mode*phis)
    return HZ



"""
Interpolate scalars on 2D poloidal plane

returns:
    values: interpolated values, values[var, harmonic, element, is, it]
"""
def interp_scalars(values, vertex, size, n_sub):
    # Multiply values[var,order,harm,vertex,element] with
    # size[order, vertex, element] and bf[order, vertex, s, t]
    return np.einsum('lihjk,ijk,ijmn->lhkmn',
                        values[:,:,:,vertex],
                        size, bf(n_sub))


"""
Interpolate scalars on 2D planes * n_planes
"""
def interp_scalars_3D(values, vertex, size, n_sub, HZ):
    values = interp_scalars(values, vertex, size, n_sub)
    return np.einsum('lhkmn,hp->lpkmn', values, HZ)


"""
Create a grid of nsub**2 points per element, at phis positions
return points and connectivity matrix
"""
def create_grid(x, vertex, size, n_elements, n_sub=4, phis=[0], n_plane=1, periodic=False):
    RZ     = grid_2D(x, vertex, size, n_sub)

    # Create connectivity data
    # Calculate 2D connectivity first
    # For each element, calculate the number of the lowest point
    # Create (n_sub-1)**2 quadrangles
    n_points = n_elements*(n_sub**2) # number of points in one plane
    n_cells  = n_elements*((n_sub-1)**2) # Number of cells in one plane
    if (n_plane > 1): # Create a volume
        if (periodic):
            n_cells_tor = n_plane
        else:
            n_cells_tor = n_plane - 1
    else:
        n_cells_tor = 1

    xyz = np.zeros((n_points*n_plane,3))
    for i in range(n_plane):
        xyz[i*n_points:(i+1)*n_points,:] = np.reshape(
            np.stack((RZ[:,:,:,0]*np.cos(phis[i]),
                      RZ[:,:,:,0]*np.sin(phis[i]),
                      RZ[:,:,:,1]), axis=-1),
            (-1,3))

    # The base block in a 2D plane
    block = np.zeros((n_sub-1,n_sub-1,4), dtype=np.int32)
    for j in range(n_sub-1):
        for k in range(n_sub-1):
            block[j,k,:] = [n_sub*j    +k  ,n_sub*(j+1)+k,
                            n_sub*(j+1)+k+1,n_sub*j    +k+1]

    i_start = np.arange(0,n_points, n_sub**2, dtype=np.int32)
    ien = np.reshape(i_start[:,np.newaxis,np.newaxis,np.newaxis]+
                     block[np.newaxis,:,:,:], (-1,4))

    # Define only _within_ an element for now
    if (n_plane > 1):
        ien_2D = ien
        ien = np.zeros((n_cells*n_cells_tor,9), dtype=np.int32)
        ien[:,0] = 8
        ien[:,1:9] = np.tile(ien_2D, (n_cells_tor,2))
        for i in range(len(phis)-1):
            ien[i*n_cells:(i+1)*n_cells,1:9] += np.concatenate(([n_points*i]*4,[n_points*(i+1)]*4))
        if (periodic):
            i = len(phis)
            ien[i*n_cells:(i+1)*n_cells,1:9] += np.concatenate(([n_points*i]*4,[0]*4))

        n_cells = n_cells * n_cells_tor
    else:
        ien = np.insert(ien, 0, 4, axis=1)

    return (xyz, ien)


"""
Calculate values of the basis functions at positions s and t
Optionally put many values of s and t at once as numpy arrays.
Dimension 0: order
Dimension 1: vertex
optional dimension 2, 3: position s, t
"""
def basis_functions(s,t):
    return np.asarray([
        [ (-1 + s)**2*(1 + 2*s)*(-1 + t)**2*(1 + 2*t),
         -(s**2*(-3 + 2*s)*(-1 + t)**2*(1 + 2*t)),
          s**2*(-3 + 2*s)*t**2*(-3 + 2*t),
         -((-1 + s)**2*(1 + 2*s)*t**2*(-3 + 2*t))],
        [ 3*(-1 + s)**2*s*(-1 + t)**2*(1 + 2*t),
         -3*(-1 + s)*s**2*(-1 + t)**2*(1 + 2*t),
         3*(-1 + s)*s**2*t**2*(-3 + 2*t),
         -3*(-1 + s)**2*s*t**2*(-3 + 2*t)],
        [ 3*(-1 + s)**2*(1 + 2*s)*(-1 + t)**2*t,
         -3*s**2*(-3 + 2*s)*(-1 + t)**2*t,
          3*s**2*(-3 + 2*s)*(-1 + t)*t**2,
         -3*(-1 + s)**2*(1 + 2*s)*(-1 + t)*t**2],
        [ 9*(-1 + s)**2*s*(-1 + t)**2*t,
         -9*(-1 + s)*s**2*(-1 + t)**2*t,
          9*(-1 + s)*s**2*(-1 + t)*t**2,
         -9*(-1 + s)**2*s*(-1 + t)*t**2]])


"""
Calculate basis functions at n_sub**2 points
"""
def bf(n_sub):
    # Get the basis functions at each of the points
    lin = np.linspace(0.0, 1.0, n_sub, dtype=prec)
    s  = np.tensordot(lin, [1]*n_sub, axes=0)
    t  = s.transpose()
    return basis_functions(s, t)
