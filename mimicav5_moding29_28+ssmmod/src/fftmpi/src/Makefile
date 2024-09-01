# Makefile for fftMPI = 2d and 3d parallel FFTs
# type "make help" for options

SHELL = /bin/sh

# ----------------------------------------
# should only need to change this section
# compiler/linker settings
# ----------------------------------------

CC =		mpicxx
CCFLAGS =	-g -O3 
SHFLAGS =	-fPIC
ARCHIVE =       ar
ARCHFLAGS =     -rc
SHLIBFLAGS =	-shared

# files

LIB2D =   libfft2dmpi.a
SHLIB2D = libfft2dmpi.so
SRC2D =	fft2d.cpp remap2d.cpp fft2d_wrap.cpp remap2d_wrap.cpp \
	memory.cpp error.cpp
INC2D =	ffttype.h fft2d.h remap2d.h pack2d.h fft2d_wrap.h remap2d_wrap.h \
	memory.h error.h
SRC90 = fft2d_wrap.f90

# memory, 1d FFT, precision settings

MEM_INC =	-DFFT_MEMALIGN=64

tbb = 		no
TBB =		$(shell echo $(tbb) | tr a-z A-Z)

ifeq ($(TBB),NO)
  TBB_INC = -DFFT_INTEL_NO_TBB
endif

fft = 		kiss
FFT =		$(shell echo $(fft) | tr a-z A-Z)
FFT_INC =       -DFFT_$(FFT)

ifeq ($(FFT),MKL)
  FFT_INC = -DFFT_$(FFT) -I$(MKL_INCLUDE)
endif

p = 		DOUBLE
P =       	$(shell echo $(p) | tr a-z A-Z)

ifeq ($(P),SINGLE)
  PRECISION = -DFFT_SINGLE
else
  PRECISION = -DFFT_DOUBLE
endif

EXTRA_INC = $(MEM_INC) $(TBB_INC) $(FFT_INC) $(PRECISION)

# 1d FFT libraries used for shared fftMPI lib build

ifeq ($(FFT),FFTW2)
  FFT1D = -lfftw
else ifeq ($(FFT),FFTW3)
  ifeq ($p,single)
    FFT1D = -lfftw3f
  else
    FFT1D = -lfftw3
  endif
else ifeq ($(FFT),FFTW)
  ifeq ($p,single)
    FFT1D = -lfftw3f
  else
    FFT1D = -lfftw3
  endif
else ifeq ($(FFT),MKL)
  ifeq ($c,intel)
    FFT1D = -L$(MKL_LIB) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
  else
    FFT1D = -L$(MKL_LIB) -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
  endif
else
  FFT1D =
endif

# targets

OBJSFFT = $(patsubst %.f90,$(BUILDDIR)/%.o,$(SRC90))
OBJ2D = $(SRC2D:.cpp=.o)

all-fft: $(SHLIB2D) $(LIB2D) $(OBJSFFT)

#$(SHLIB2D):
#	$(MAKE) $(SHARED2D) fft=$(fft) p=$(p)

#lib2d:
#	$(MAKE) static2d fft=$(fft) p=$(p)

$(SHLIB2D):	$(OBJ2D)
		$(CC) $(CCFLAGS) $(SHFLAGS) $(SHLIBFLAGS) -o $(SHLIB2D) $(OBJ2D) $(FFT1D) 

$(LIB2D):	$(OBJ2D)
		$(ARCHIVE) $(ARCHFLAGS) $(LIB2D) $(OBJ2D)

clean-fft: 
	rm -f *.o *.pyc

$(OBJ2D) : %.o :%.cpp
	$(CC) $(CCFLAGS) $(SHFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@)

$(OBJSFFT) : $(BUILDDIR)/%.o:%.f90
	$(FC) $(FFLAGS) $(EXTRA_INC) -c $< -o $(BUILDDIR)/$(notdir $@)
