
export NETCDF='/software/meteo/xenial/x86_64/netcdf/4.4.1.1-gcc'
NETCDFDIR= $(NETCDF)

EXECNAME = cloud_tracking.exe

FC = gfortran
FFLAGS = -ffree-line-length-512 -I$(NETCDFDIR)/include 
NETCDFFLAGS = -L$(NETCDFDIR)/lib -L$(NETCDFDIR)/lib64 -lnetcdff

SRC=clouds_ct.f90 shared_ct.f90 cloud_search.f90 netcdf_interface.f90 cloud_tracking_sa.f90

OBJSF90=$(patsubst %.f90,%.o,$(SRC))

all: $(OBJSF90)
	$(FC) -o $(EXECNAME) ./*.o $(NETCDFFLAGS)

$(OBJSF90) : %.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $(notdir $@)
	
clean:
	rm $(EXECNAME) *.o *.mod
