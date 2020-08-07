# Contents of this file are included in palm's makefile if the simple install
# process (palm_simple_install) is used.
# Please adjust compiling parameters as required for your system.

PROG     = palm
F90      = mpif90

COPT     = -cpp -DMPI_REAL=MPI_DOUBLE_PRECISION -DMPI_2REAL=MPI_2DOUBLE_PRECISION -D__parallel  -D__lc  -D__netcdf4_parallel  -D__netcdf -D__netcdf4
F90FLAGS = -xCORE-AVX512 -fma -O3  -no-prec-sqrt -ip -convert little_endian  -unroll-aggressive -no-prec-div -fp-model fast=2   -I$(PARALLEL_NETCDF_INSTALL_ROOT)/include -I$(NETCDF_FORTRAN_INSTALL_ROOT)/include
LDFLAGS  = -xCORE-AVX512 -fma -O3  -no-prec-sqrt -ip -convert little_endian  -unroll-aggressive -no-prec-div -fp-model fast=2   -L$(PARALLEL_NETCDF_INSTALL_ROOT)/lib  -lpnetcdf -L$(NETCDF_FORTRAN_INSTALL_ROOT)/lib  -lnetcdff

#
# Currently Loaded Modules:
#  1) intel/19.0.4   2) hpcx-mpi/2.4.0   3) intel-mkl/2019.0.4   4) StdEnv   5) parallel-netcdf/1.8.0   6) hdf5/1.10.4-mpi   7) netcdf-fortran/4.4.4



# -ftz -xCORE-AVX512 -no-prec-div -no-prec-sqrt -ip -convert little_endian

