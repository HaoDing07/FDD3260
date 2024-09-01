
export NETCDF='/software/opt/bionic/x86_64/netcdf/4.7.0-gcc'
NETCDFDIR= $(NETCDF)

EXECNAME = slices_x.exe

FC = gfortran
FFLAGS = -ffree-line-length-1024 -I$(NETCDFDIR)/include -fdefault-real-8 
FFLAGS := $(FFLAGS) -g -O0 -fbounds-check -ffpe-trap=invalid,overflow,zero -fbacktrace 
NETCDFFLAGS = -L$(NETCDFDIR)/lib -L$(NETCDFDIR)/lib64 -lnetcdff

SRC=shared_pr.f90 netcdf_interface.f90 diagnostics.f90 slices.f90

OBJSF90=$(patsubst %.f90,%.o,$(SRC))

all: $(OBJSF90)
	$(FC) -o $(EXECNAME) ./*.o $(NETCDFFLAGS)

$(OBJSF90) : %.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $(notdir $@)
	
clean:
	rm $(EXECNAME) *.o *.mod
