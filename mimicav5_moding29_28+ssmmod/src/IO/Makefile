SRC=netcdf_define.f90 netcdf_slice.f90 netcdf_interface.f90 check.f90 diagnostics.f90 initial.f90 output.f90 outputpack.f90 restarting.f90 setrun.f90

OBJSF90 = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRC))

all-io: $(OBJSF90)

clean-io: 	
	rm -f ./*.mod

$(OBJSF90) : $(BUILDDIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(NETCDFFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@) 
