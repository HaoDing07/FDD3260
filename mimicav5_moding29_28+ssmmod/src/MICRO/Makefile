SRC=micropack.f90 activation.f90 cnt.f90 freezing.f90 kessler.f90 sbpack.f90 sb_warm.f90 sb_ice.f90 precipitation.f90 micro_diag.f90  micro_interface.f90  

OBJSF90 = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRC))

all-micro: $(OBJSF90)

clean-micro: 
	rm -f ./*.mod

$(OBJSF90) : $(BUILDDIR)/%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@)
