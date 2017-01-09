#!/bin/bash

module purge
module load intel/pe-xe-2017--binary \
          intelmpi/2017--binary \
          mkl/2017--binary \
          fftw/3.3.4--intelmpi--2017--binary \
          zlib/1.2.8--gnu--6.1.0 \
          szip/2.1--gnu--6.1.0 \
          hdf5/1.8.17--intel--pe-xe-2017--binary \
          lapack/3.6.1--intel--pe-xe-2017--binary \
          blas/3.6.0--intel--pe-xe-2017--binary 

export PASTIX_LIBDIR=/marconi/home/userexternal/glatu000/jorek_tools/pastix_3184_knl_17/install
export MUMPS_HOME=/marconi/home/userexternal/glatu000/jorek_tools/mumps-5.0.2_knl
export LANG=C
export JOREK_HOST=marconi
export compilethreads=32
export MPIRUN="mpirun -genv I_MPI_DEBUG=4 -genv KMP_AFFINITY=verbose,scatter --map-by ppr:1:socket:pe=$OMP_NUM_THREADS -np"
export BATCHCOMMAND="qsub"
export http_proxy=http://proxy:3128 
