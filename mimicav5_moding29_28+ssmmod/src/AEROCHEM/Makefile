SRC=aeroactivate.f90 aeroscavenge.f90 aerosol.f90 aq.f90 chemfront.f90 chemscavenge.f90 chemtrop.f90 nucleate.f90 pdadv.f90 solidtrop.f90

OBJSF90 = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRC))

all-aero: $(OBJSF90)

clean-aero: 
	rm -f ./*.mod

$(OBJSF90) : $(BUILDDIR)/%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@)
