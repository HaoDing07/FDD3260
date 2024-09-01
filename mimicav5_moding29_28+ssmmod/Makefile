##
## === Define files and directories
##

SRCDIR  = $(MIMICA)/src
INCDIR  = $(MIMICA)/INCLUDE
FFTMPIDIR = $(MIMICA)/src/fftmpi
BUILDDIR= $(PWD)/build
NETCDFDIR= $(NETCDF_FORTRAN_BASE)
EXECNAME = mimicav5_$(VERSION).exe

SRCDIR_TESTS = $(MIMICA)/pFUnitTests
BUILDDIR_TESTS = $(PWD)/pFUnitTests 

##
## === Compiler specific options (Intel)
##

ifeq ($(CMPLER),INTEL)

FC = ifort
ifeq ($(SPMD),TRUE)
FC = mpif90
endif

FFLAGS = -I$(INCDIR) -I$(BUILDDIR) $(LIBFLAGS) -cpp -r8 -lstdc++ -module $(BUILDDIR)

ifeq ($(DEBUG),TRUE)
FFLAGS := $(FFLAGS) -g -O1 -fpe0 -ftrapuv -fp-model strict -fp-stack-check -CB -CU -check nooutput_conversion -debug all -traceback
else
## FFLAGS := $(FFLAGS) -O3 -fpe3 -standard-semantics -assume nostd_mod_proc_name -no-prec-div -xHost -nowarn -heap-arrays
FFLAGS := $(FFLAGS) -O3 -fpe3 -ftz -ip -xHost -no-prec-div -no-prec-sqrt -nowarn -heap-arrays
endif
ifeq ($(PROF),TRUE)
FFLAGS := $(FFLAGS) -O0 -g -pg
LIBFLAGS := $(LIBFLAGS) -pg
endif
ifeq ($(SPMD),TRUE)
FFLAGS := $(FFLAGS) -shared-intel -mcmodel=medium -lmpi
endif

else

##
## === Compiler specific options (Gfortran - default)
##

FC = gfortran
ifeq ($(SPMD),TRUE)
FC = mpif90
endif

FFLAGS = -I$(INCDIR) -J$(BUILDDIR) $(LIBFLAGS) -cpp -ffree-line-length-none -fdefault-real-8

ifeq ($(DEBUG),TRUE)
FFLAGS := $(FFLAGS) -O0 -g -finit-real=snan -fimplicit-none -ffpe-trap=invalid,overflow,zero -fbounds-check -fbacktrace 
else
FFLAGS := $(FFLAGS) -mtune=native -O3 -ffast-math -ftree-vectorize
endif
ifeq ($(PROF),TRUE)
FFLAGS := $(FFLAGS) -O0 -g -pg
LIBFLAGS := $(LIBFLAGS) -pg
endif
ifeq ($(SPMD),TRUE)
FFLAGS := $(FFLAGS) -mcmodel=medium -lmpi -fmax-stack-var-size=2048
endif

endif

##
## === Generic external libraries
##

LIBFLAGS := $(LIBFLAGS) -lm
NETCDFFLAGS = -L$(NETCDFDIR)/lib -lnetcdff -I$(NETCDFDIR)/include

FFTMPIFLAGS = -I$(FFTMPIDIR)/src -L$(FFTMPIDIR)/src -lfft2dmpi
FFTFLAGS = $(FFTW_LIB) $(FFTW_MPI_LIB) $(FFTW_INC)
ifneq ($(CMPLER),INTEL)
FFTFLAGS := $(FFTFLAGS) -std=f2003
endif

export FC
export FFLAGS
export FFTFLAGS
export FFTMPIFLAGS
export NETCDFFLAGS
export LIBFLAGS
export BUILDDIR
export SRCDIR_TESTS
export BUILDDIR_TESTS

##
## === rules
##
     
all:
	mkdir -p $(BUILDDIR)
	cd $(SRCDIR) ; $(MAKE)
	$(FC) -o $(EXECNAME) $(BUILDDIR)/*.o $(NETCDFFLAGS) $(FFTFLAGS) $(FFTMPIFLAGS) $(LIBFLAGS)

clean:
	rm -rf $(BUILDDIR)/* ; rm -f mimicav5_*.exe
	cd $(SRCDIR) ; $(MAKE) clean
