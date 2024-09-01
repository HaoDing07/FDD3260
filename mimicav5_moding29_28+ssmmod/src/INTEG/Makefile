SRC=integpack.f90 windadv.f90 momentum.f90 scalars.f90 tracers.f90 prognostic_pressure.f90 time_step.f90 advcm.f90 column_step.f90 column.f90 piggy.f90

OBJSF90 = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRC))

all-integ: $(OBJSF90)

clean-integ: 
	rm -f ./*.mod

$(OBJSF90) : $(BUILDDIR)/%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@)
