SRC=system_tools.f90 alloc.f90 bc.f90 advq_ppm.f90 advq_weno.f90 advq_muscl.f90 advq_quick.f90 advq_quick_nolim.f90 advq_lw.f90 advq_lw_nolim.f90 gradients.f90 averages.f90 thermodynamics.f90 lagrange_interp.f90 linear_interp.f90 randpack.f90 define_grid.f90 hygroscopic_growth.f90 surface_var.f90 timevar.f90 vof.f90 nesting.f90 funcpack.f90 advection.f90

SRC_FFT=mpi_prepare.f90

OBJSF90 = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRC))

OBJSFFT = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRC_FFT))

all-tools: $(OBJSF90) $(OBJSFFT)

clean-tools: 
	rm -f ./*.mod

$(OBJSF90) : $(BUILDDIR)/%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@)
        
$(OBJSFFT) : $(BUILDDIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(FFTFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@)
