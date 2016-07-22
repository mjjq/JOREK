title: Installation

## Dependencies
* A working JOREK installation
* HDF5, set paths and `USE_HDF5 = 1` in `Makefile.inc`
* MUMPS, set paths and `USE_MUMPS = 1` in `Makefile.inc`

### HDF5
Can probably be obtained from your module system or distribution package manager.
Parallel HDF5 is needed in JOREK2_particles, so it might be needed to compile this separately.
See [INSTALL_parallel](https://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/release_docs/INSTALL_parallel) for more details.
It is important to also set the `--enable-fortran` option.

### MUMPS
* Download the latest version from [the MUMPS page](http://mumps.enseeiht.fr/) or download [MUMPS_5.0.2.tar.gz](http://mumps.enseeiht.fr/MUMPS_5.0.2.tar.gz)
* Copy the relevant file from the `Make.inc` directory to `Makefile.inc`
* If you want to compile with (PT)SCOTCH, follow the instructions on [the JOREK wiki](http://jorek.eu/wiki/doku.php?id=compiling#scotch).
    * Do not download the version with `libesmumps`.
* Set the paths in your `Makefile.inc` and select which one of PORD, SCOTCH or METIS you want to use.
* `make alllib` compiles `c z s d` versions. You need at least `s` and `d`.
* See the `INSTALL` file for more info.


## Compilation
You can compile single-threaded
```bash
make jorek2_particles
```
or use multiple threads
```bash
make -j8 jorek2_particles
```
In this case you might need to restart the compilation once, as dependency information is not yet fully correct.

### Diagnostics
Compile some of the diagnostics
```bash
make count_particles_vtk
make dump_particles_vtk
make particle_flux_coordinates
make particle_flux_coordinate_diffusion
make project_particles_vtk
```
