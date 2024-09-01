SRC = modelctl.f90 mimica.f90

OBJSF90 = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRC))

all-main: $(OBJSF90) 

clean-main: 
	rm -f ./*.mod

$(OBJSF90) : $(BUILDDIR)/%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@)
