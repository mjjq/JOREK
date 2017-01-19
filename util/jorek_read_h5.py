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

import h5py
import numpy as np


"""
Class encapsulating some read data from a restart file

input arguments:
    Variables: a dict of {key: name} for all variables you want
    i_plane: which plane to interpolate at. If -1 make a 3D field.
"""
class fields:
    var_names = ["psi", "u", "j", "w", "rho", "T", "v_par"] # rest is ambiguous
    def read(self, filename, variables=''):
        if (isinstance(variables, str)):
            self.variables = self.variablesToDict(variables)
        else:
            self.variables = variables
        self.var_nums = list(self.variables.keys())

        with h5py.File(filename, 'r') as hf:
            self.n_var        = hf.get('n_var')[0]
            self.n_period     = hf.get('n_period')[0]
            self.n_tor        = hf.get('n_tor')[0]
            self.n_vertex_max = hf.get('n_vertex_max')[0]
            self.n_elements   = hf.get('n_elements')[0]
            self.vertex       = np.array(hf.get('vertex'), dtype=np.int32)-1
            self.x            = np.array(hf.get('x'))
            self.size         = np.array(hf.get('size'))
            if (len(self.var_nums) > 0):
                self.values       = np.array(hf.get('values'))[self.var_nums,:,:]

    """
    Convert a list of variables like 1,2,3
    to the dict we expect
    """
    def variablesToDict(self, var_string):
        dict = {}
        for var in var_string.split(','):
            if (var.strip().isdigit()):
                i = int(var.strip())
                dict[i] = self.var_names[i-1]
        return dict

    """
    Create a set of points and cells and interpolate values 
    """
    def interp_to_mesh(self, n_sub=4, n_plane=40, phi=[0,2*np.pi], without_n0_mode=False):
        if not hasattr(self, 'n_tor'):
            print("ERROR: no file read yet")
            return
        # Get the basis functions at each of the points
        lin = np.linspace(0.0, 1.0, n_sub, dtype=np.double)
        s  = np.tensordot(lin, [1]*n_sub, axes=0)
        t  = s.transpose()
        bf = basis_functions(s, t)

        # If we do not have a range in phi make only one plane
        if (isinstance(phi,int)):
            n_plane = 1
            phis = phi
        elif (n_plane == 1):
            phis = np.asarray([phi[0]])
        else:
            phis = np.linspace(phi[0],phi[1],num=n_plane,endpoint=False)

        # Setup toroidal coefficients for each plane and toroidal harmonic
        HZ = np.zeros((self.n_tor,n_plane))
        for i in range(self.n_tor):
            mode = np.floor((i+1)/2)*self.n_period
            if (i == 0):
                if (not without_n0_mode):
                    HZ[i,:] = 1
            elif (i % 2 == 0):
                HZ[i,:] = np.sin(mode*phis)
            elif (i % 2 == 1):
                HZ[i,:] = np.cos(mode*phis)
        
        # Calculate RZ for all of the elements (dimension 0) for each of the s
        # positions (dimension 1) for each of the t positions (dimension 2)
        # Multiply x[var, order, vertex, element] on the last two dimensions
        # with size[order, vertex, element]*bf[order, vertex, s, t]
        # See http://stackoverflow.com/questions/26089893/understanding-numpys-einsum
        x = np.einsum('lijk,ijk,ijmn->kmnl', self.x[:,:,self.vertex], self.size, bf)
        # x[element,is,it,var]

        # Multiply values[var,order,harm,vertex,element] with
        # size[order, vertex, element] and bf[order, vertex, s, t]
        if (len(self.var_nums) > 0):
            scalars = np.einsum('lihjk,ijk,ijmn->lhkmn',
                                self.values[:,:,:,self.vertex],
                                self.size, bf)

        # TODO: use numba for GPU / vectorization?
        #print(scalars)
        return (x, scalars)


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

