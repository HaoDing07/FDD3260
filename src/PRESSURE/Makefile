SRC_F90_FFT=fft_tools.f90 fft_interface.f90 fft_interface_open.f90 fft_interface_nnf.f90 fftmpi_interface.f90 pressure_solver.f90 

OBJSFFT = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRC_F90_FFT))

all-pressure: $(OBJSFFT)

clean-pressure: 
	rm -f ./*.mod

$(OBJSFFT) : $(BUILDDIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(FFTWFLAGS) $(FFTMPIFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@)
