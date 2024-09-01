SRC=util.f90 mo_opac_optical_properties.f90 mo_aerorad_lut.f90 lagrange.f90 nudge.f90 rad_corkds.f90 rad_cldwtr.f90 rad_rndnmb.f90 rad_solver.f90 scalar_sources.f90 surfacepack.f90 subgrid.f90 getcape_Bryan.f90 column_phys.f90
SRCLAST=rad_fuliou.f90 radiag.f90

OBJSF90 = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRC) $(SRCLAST))

all-phys: $(OBJSF90)

clean-phys: 
	rm -f ./*.mod

$(OBJSF90) : $(BUILDDIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(NETCDFFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@)
