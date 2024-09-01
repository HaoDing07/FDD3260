SRCF=gridno.f90
SRCFIRST=shared_data.f90 radia.f90
SRCTYPE=typedef_aerosol_2.f90 typedef_aq.f90 typedef_gas.f90 typedef_hydrometeor.f90 typedef_lagrange.f90 typedef_nuclei.f90 typedef_rad.f90 typedef_solid.f90 typedef_state.f90 typedef_surf.f90 typedef_svalue.f90 typedef_thermo.f90 typedef_turbu.f90 typedef_wind.f90 typedef_tendencies.f90 typedef_pressure.f90 typedef_diag.f90
SRCSHARED=shared_aq.f90 shared_gas.f90 shared_hydro.f90 shared_lagrange.f90 shared_rad.f90 shared_solid.f90 shared_state.f90 shared_surf.f90 shared_aerosol_2.f90 shared_nuclei.f90 shared_thermo.f90 shared_turbu.f90 shared_wind.f90 shared_tendencies.f90 shared_pressure.f90 shared_diag.f90 shared_all.f90
SRCLAST=mpicomm.f90

OBJSF = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRCF))

OBJSF90 = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRCFIRST) $(SRCTYPE) $(SRCSHARED) $(SRCLAST))

all-types: $(OBJSF) $(OBJSF90)

clean-types: 
	rm -f ./*.mod ./gridno.f90

$(OBJSF90) : $(BUILDDIR)/%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@)

$(OBJSF) : $(BUILDDIR)/%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@)

$(BUILDDIR)/%.mod: $(BUILDDIR)/%.o
	@true
