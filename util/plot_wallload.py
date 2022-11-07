#!/usr/bin/env python
"""
Script for postprocessing and plotting the wall loads evaluated with the particle tracer.

See particles/mod_wall_collision.f90 on how the wall loads are evaluated. This script is used
to plot and postprocess the results. To run this script, you'll need the wall input (HDF5
file containing the wall triangles) and either/both the resulting particle restart file
(where the wall IDs are stored in i_elm field, see mod_wall_collision.f90) and
the wall load file that can be exported with a routine in mod_wall_collision.f90. Place those
files in the folder together with this script, type in the filenames below and just run this
script.

This script needs pyvista to work, but you can comment wherever it is used to plot
some rudimentary results. However, pyvista is required to convert the data to a vtk file
(this is done by default and stored in wallload.vtk).

The 3D plot is by default shown in matplotlib, but you can view the wall loads interactively
by setting showpyvista = True. All in all, you should modify this file to suit your exact
needs.
"""
import numpy as np
import h5py
import pyvista as pv

## Settings ##

# Make plots or just output the 0D quantities
plot2d = True
plot3d = True

# Axis (R,z) coordinates (these determine poloidal angle in plots)
r0 = 6.2
z0 = 0.1

# Wall input data
wallin   = "iterwall_offset20cm.h5"

# Wall loads
wallload = "wallload.h5"

# Particle restart file containing the output (or None)
partout  = "part_out.h5"
partout = None

# Group to be plotted
group = "001"

# Define a contour along the wall for plotting losses along it
limiter = np.array([[4,4,4.2,4.9,5.7,6.5,6.8],
                    [3,4,4.2,4.6,4.4,3.8,3.5]]).T

# Choose whether to show interactive view using pyvista or final product with pyplot
showpyvista = False

## End of settings ##

# Hack to convert e -> x when displaying numbers
def as_si(x, ndp):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))


with h5py.File(wallin,'r') as h5:
    data = h5['nodes'][:]
    ntriangle = h5['ntriangle'][0]

    wr = np.sqrt(data[0::3]**2 + data[1::3]**2)
    wz = data[2::3]

    verts = []
    faces = []
    wd    = np.zeros((ntriangle,))
    wr    = np.zeros((ntriangle,))
    wz    = np.zeros((ntriangle,))
    wphi  = np.zeros((ntriangle,))
    wpol  = np.zeros((ntriangle,))
    warea = np.zeros((ntriangle,))
    for i in range(ntriangle):

        # Compute area for each node
        a = data[i*9+3:i*9+6] - data[i*9+0:i*9+3] 
        b = data[i*9+6:i*9+9] - data[i*9+0:i*9+3]
        c = np.cross(a,b)
        warea[i] = np.sqrt( np.sum(c*c) ) / 2

        # Compute (r,z,tor,pol) coordinates for each element center (for plotting)
        x = np.sum(data[i*9+0:i*9+7:3]) / 3
        y = np.sum(data[i*9+1:i*9+8:3]) / 3
        z = np.sum(data[i*9+2:i*9+9:3]) / 3

        wr[i]   = np.sqrt(x**2 + y**2)
        wz[i]   = z
        wphi[i] = np.mod(np.arctan2( y, x ), 2*np.pi)
        wpol[i] = np.arctan2( z - z0, wr[i] - r0 )

        # Compute distance along the contour
        dist = 100
        idx  = 0
        for j in range(limiter.shape[0]-1):
            dist0 = ( np.abs( (limiter[j+1,0]-limiter[j,0])*(limiter[j,1]-wz[i]) - (limiter[j,0]-wr[i])*(limiter[j+1,1]-limiter[j,1]) )
                      / np.sqrt( (limiter[j+1,0]-limiter[j,0])**2 + (limiter[j+1,1]-limiter[j,1])**2 ) )
            if dist0 < dist:
                dist = dist0
                idx  = j

        for j in range(idx):
            wd[i] += np.sqrt( (limiter[j+1,0]-limiter[j,0])**2 + (limiter[j+1,1]-limiter[j,1])**2 )

        j = idx
        wd[i] += ( (wr[i] - limiter[j,0]) * (limiter[j+1,0]-limiter[j,0]) + (wz[i] - limiter[j,1]) * (limiter[j+1,1]-limiter[j,1]) ) \
                 / np.sqrt( (limiter[j+1,0]-limiter[j,0])**2 + (limiter[j+1,1]-limiter[j,1])**2 )

        # Collect vertices and faces for pyvista
        verts += [data[i*9+0],data[i*9+1],data[i*9+2]]
        verts += [data[i*9+3],data[i*9+4],data[i*9+5]]
        verts += [data[i*9+6],data[i*9+7],data[i*9+8]]
        faces += [[3, i*3 + 0, i*3 + 1, i*3 + 2]]

    wallmesh = pv.PolyData(verts, faces)
    wallmesh.cell_data['pload'] = np.zeros((ntriangle,))
    wallmesh.cell_data['eload'] = np.zeros((ntriangle,))

    del data, verts, faces

with h5py.File(wallload,'r') as h5:
    wallid   = h5[group]["wallid"][:] - 1 # Fortran indexing starts at 1, python at 0
    enedepot = h5[group]["energydepot"][:]
    prtdepot = h5[group]["particledepot"][:]

    wetted_area = np.sum(warea[wallid])
    peak_pload  = np.amax(prtdepot/warea[wallid])
    peak_eload  = np.amax(enedepot/warea[wallid])

    wallmesh.cell_data['pload'][wallid] = prtdepot/warea[wallid]
    wallmesh.cell_data['eload'][wallid] = enedepot/warea[wallid]
    wallmesh.save('wallload.vtk')
    
    print("Wetted area: "        + "{:.2e}".format(wetted_area) + " m^2")
    print("Peak particle load: " + "{:.2e}".format(peak_pload)  + " prt/m^2")
    print("Peak energy load: "   + "{:.2e}".format(peak_eload)  + " J/m^2")

if partout is not None:
    with h5py.File(partout,'r') as h5:
        wallid_prt = h5["groups/"+group+"/i_elm"][:]
        r_prt   = h5["groups/"+group+"/x"][:,0]
        z_prt   = h5["groups/"+group+"/x"][:,1]
        phi_prt = np.mod(h5["groups/"+group+"/x"][:,2], 2*np.pi)
        pol_prt = np.arctan2( z_prt - z0, r_prt - r0 )
        weight  = h5["groups/"+group+"/weight"][:]

        lost = wallid_prt < 0
        wallid_prt = -wallid_prt - 1

        wettedid = np.unique(wallid_prt[lost])
        wetted_area_prt = np.sum(warea[wettedid])

        prtdepot_prt = np.zeros((wettedid.size,))
        for i,e in enumerate(wettedid):
            prtdepot_prt[i] = np.sum(weight[e == wallid_prt])

        peak_pload_prt = np.amax(prtdepot_prt / warea[wettedid])

        print("\nBased on particle restart file:")
        print("Wetted area: "        + "{:.2e}".format(wetted_area_prt) + " m^2")
        print("Peak particle load: " + "{:.2e}".format(peak_pload_prt)  + " prt/m^2")


if plot2d or plot3d:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    cm = 1/2.54
    params = {'legend.fontsize': 8,
              'axes.labelsize':  8,
              'axes.titlesize':  8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'font.size' : 8,
              'text.usetex' : True}
    plt.rcParams.update(params)


if plot2d:
    import matplotlib


    fig1 = plt.figure()
    gs = GridSpec(1,2,figure=fig1,width_ratios=[1, 2], right=0.85)
    s1 = fig1.add_subplot(gs[0])
    s2 = fig1.add_subplot(gs[1])

    # In case pyvista is not available
    pyvista = True
    if not pyvista:
        s1.scatter( wr, wz, 1, 'black' )

    else:
        planemesh   = pv.Plane(center=(r0,0,z0),direction=(0,1,0),i_size=r0*2,j_size=10).triangulate()
        wallcontour = wallmesh.intersection(planemesh)[0]
        pts = wallcontour.points
        idx = wallcontour.lines

        i = 0
        while i < idx.size:
            nk = idx[i]
            xyz = pts[idx[i+1:i+1+nk],:]
            s1.plot(xyz[:,0], xyz[:,2], color='black')
            i = i + nk + 1

    if partout is not None:
        h1 = s1.scatter(r_prt[lost], z_prt[lost], 1, 'red')
        h2 = s1.scatter(r_prt[~lost], z_prt[~lost], 1, 'blue')

        s1.legend([h1, h2], ["Lost", "Confined"])

        # Uncomment to plot marker positions
        #s2.scatter(phi_prt[lost], pol_prt[lost], 1, 'red')
        h3 = s2.scatter(wphi[wettedid], wpol[wettedid], 2, prtdepot_prt/warea[wettedid], cmap='Reds')

        title = r"prt/m$^2$"
        s2.annotate("Wetted area: "        + r"${0:s}$".format(as_si(wetted_area_prt,2)) + r" m$^2$",     (0.1,-1.4) )
        s2.annotate("Peak particle load: " + r"${0:s}$".format(as_si(peak_pload_prt,2))  + r" prt/m$^2$", (0.1,-1.7) )
    else:
        s1.scatter(wr[wallid], wz[wallid], 1, 'red')
        h3 = s2.scatter(wphi[wallid], wpol[wallid], 2, enedepot/warea[wallid], cmap='Reds')
        
        title = r"J/m$^2$"
        s2.annotate("Wetted area: "        + r"${0:s}$".format(as_si(wetted_area,2)) + r" m$^2$",     (0.1,-1.4) )
        s2.annotate("Peak particle load: " + r"${0:s}$".format(as_si(peak_pload,2))  + r" prt/m$^2$", (0.1,-1.7) )
        s2.annotate("Peak energy load: "   + r"${0:s}$".format(as_si(peak_eload,2))  + r" J/m$^2$",   (0.1,-2.0) )

    cax = plt.colorbar(h3, ax=s2, location='top')
    cax.set_label(title)
        
    s2.set_xlim(0, 2*np.pi)
    s2.set_ylim(-np.pi, np.pi)

    s2.set_xlabel("Toroidal angle [rad]")
    s2.set_ylabel("Poloidal angle [rad]")

    s2.set_xticks([0, np.pi, 2*np.pi])
    s2.set_xticklabels([r"0",r"$\pi$",r"$2\pi$"])

    s2.set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
    s2.set_yticklabels([r"$-pi$", r"$-pi/2$ (Divertor)", r"$0$ OMP", r"$\pi/2$ (Top)",r"$\pi$"])

    s2.yaxis.set_label_position("right")
    s2.yaxis.tick_right()

    s1.set_aspect('equal', 'box')
    s1.plot(limiter[:,0], limiter[:,1])

    
    fig2 = plt.figure()
    s1 = fig2.add_subplot(1,1,1)

    if partout is not None:
        dens,xg,yg = np.histogram2d(wphi[wettedid], wd[wettedid], bins=[90, 40], weights=prtdepot_prt/warea[wettedid])
        dens0 = np.histogram2d(wphi[wettedid], wd[wettedid], bins=[90, 40])[0]
    else:
        dens,xg,yg = np.histogram2d(wphi[wallid], wd[wallid], bins=[90, 40], weights=enedepot/warea[wallid])
        dens0 = np.histogram2d(wphi[wallid], wd[wallid], bins=[90, 40])[0]

    dens0[dens0 == 0] = 1 # Let's not divide by zero and accidentally create a black hole
    xg = xg[:-1] + (xg[1] - xg[0]) / 2
    yg = yg[:-1] + (yg[1] - yg[0]) / 2
    cmap = mpl.colormaps["Reds"].copy()
    cmap.set_bad(color=[0.9,0.9,0.9])
    h4 = s1.contourf(xg,yg,(dens/dens0).T, cmap=cmap)

    s1.set_xlabel("Toroidal angle [rad]")
    s1.set_ylabel("Distance along wall [m]")
    cax = plt.colorbar(h4, ax=s1, location='top')
    #cax.set_label("Mean load " + "$\log_{10}$" + title)
    cax.set_label("Mean load " + title)
    
    plt.tight_layout()

if plot3d:
    # Camera toroidal position
    phicam = 330 * np.pi / 180

    # Removes bottom mesh so it is not in the way of camera
    bounds = [-20, 20, -20, 20, -10, -2]
    wallmesh = wallmesh.clip_box(bounds)

    # This makes non-loaded wall elements appear grey instead red-tinted
    cmap = mpl.colormaps["Reds"].copy()
    cmap.set_bad(color=[0.9,0.9,0.9])

    if showpyvista:
        p = pv.Plotter()
    else:
        p = pv.Plotter(off_screen=True)
        p.store_image = True
        
    p.add_mesh(wallmesh, scalars='eload', cmap=cmap)
    if not showpyvista:
        p.remove_scalar_bar()
    p.camera.position = (r0*np.cos(phicam), r0*np.sin(phicam), -6)
    p.camera.focal_point = (0.9*r0*np.cos(phicam), 0.9*r0*np.sin(phicam), -0.1)
    p.show()

    if not showpyvista:
        fig3 = plt.figure(figsize=(6,6))
        s1 = fig3.add_subplot(1,1,1)
        s1.imshow(p.image)

        norm = mpl.colors.Normalize(vmin=0,vmax=peak_eload)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cax = plt.colorbar(sm,ax=s1,location='top')
        cax.set_label(r"Energy load J/m$^2$")

        s1.set_xticks([])
        s1.set_yticks([])

        plt.tight_layout()

if plot2d or plot3d:
    #fig1.savefig("wallloadscatter.png", format="png", dpi=96*2)
    #fig2.savefig("wallloadmean.png", format="png", dpi=96*2)
    #fig3.savefig("wallload3d.png", format="png", dpi=96*2)
    plt.show()
