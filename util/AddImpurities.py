#!/usr/bin/env python
# coding: utf-8

# # Import and rewrite an h5 file to add impurities as a constant or profile.

#Load packages
import numpy as np
import h5py
import sys
import select
import argparse
import os
import shutil

def main():
    parser=argparse.ArgumentParser(description='''This function adds impurities to the JOREK restart file.''')
    parser.add_argument("-c", "--constant", type=float, help="Set a constant density of impurities in JOREK units")
    parser.add_argument("-p", "--profile", type=float, help="Set the impurity profile relative to the density profile")
    parser.add_argument("-f", "--restartfile", type=str, nargs='*', help='JOREK restart file name')             
    args=parser.parse_args()                                            

    #Print a confirmation message
    print "Your input:"
    if args.constant is not None:
        print "Constant density = ", args.constant
        profile = False
    if args.profile is not None:
        print "Density profile scaled with a factor ", args.profile
        profile = True
    if (args.constant is not None and args.profile is not None):
        print "Error: Can not simultaneously set a constant level and a relative profile"
        sys.exit(0)

    #Test if input file is present
    if args.restartfile is None:
        parser.print_help()
        sys.exit(0)
    
    filename = args.restartfile[0]

    n0 = 1e20
    type_out = None
    with h5py.File(filename, 'r+') as f:
        #f.visit(print)
        
        # Obtain the values
        ds = f['values']
        offset = ds.id.get_offset()
        #print(offset)
        if (ds.chunks is None and ds.compression is None and offset > 0):
            dtype = ds.dtype
            shape = ds.shape
            arr = np.memmap(filename, mode='r+', shape=shape, offset=offset, dtype=dtype)
            if (type_out is not None and type_out is not dtype):
                print(arr.astype(type_out))
        else:
            if (type_out is None):
                arr = np.array(f.get(h5path))
            else:
                arr = np.array(f.get(h5path), dtype=type_out)
        
        # Obtain central_density
        ds = f['central_density']
        offset = ds.id.get_offset()
        #print(offset)
        if (ds.chunks is None and ds.compression is None and offset > 0):
            dtype = ds.dtype
            shape = ds.shape
            rho_0 = np.memmap(filename, mode='r', shape=shape, offset=offset, dtype=dtype)
            if (type_out is not None and type_out is not dtype):
                print(rho_0.astype(type_out))
        else:
            if (type_out is None):
                rho_0 = np.array(f.get(h5path))
            else:
                rho_0 = np.array(f.get(h5path), dtype=type_out)
        
    rho_0=float(rho_0)
    #if (float(rho_imp*rho_0*n0) != float(rho_imp_SI)) :
     #   print('Impurities not equal, check!!')
      #  print(float(rho_imp*rho_0*n0),float(rho_imp_SI))
       # print(rho_imp*rho_0*n0-rho_imp_SI)
        
    # Now array contains values(Node number, Toroidal harmonic (first one for axisymmetry), [value, s derivative, t derivative, st derivative], which variable # (1-8))
    # In python [Variable number, value_or_derivative, n_tor, node]
    # 7=rho_imp, 4=rho
    
    print "Output:"
    if profile == True : #Follow electron density profile
        arr[7,:,:,:] = arr[4,:,:,:]*args.profile
        print "Maximum impurity density = ", np.amax(arr[7,:,:,:]) * n0 * rho_0, " [m^-3]"

    elif profile == False : #Constant density
        arr[7,:,:,:] = 0.0
        arr[7,0,0,:] = args.constant
        print "Constant impurity density = ", np.amax(arr[7,0,0,:]) * n0 * rho_0, " [m^-3]"

    #After you are done
    f.close()

if __name__ == "__main__":
    main()



