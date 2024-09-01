# README #

This README file documents all the steps necessary to get your application up and running.

### What is this repository for? ###

* MIMICA is a 3D high-resolution model designed for atmospheric applications and written entirely in fortran. 
The code has been originally developed by Chien Wang at MIT, and later updated and transformed into a true 
LES model by Julien Savre at Stockholm University. Today, MIMICA is developed under the git version control 
system and is hosted on the Bitbucket web-page. For more information and technical details, please refer to 
MIMICA's manual in the ./MANUAL folder.

* The current working version of MIMICA v5 (master) has been issued on March 11th, 2019.

### How do I get set up? ###

* 1. (a) Compilation: External libraries

		MIMICA V5 is most efficiently run using the FFT based pressure solver. This latter uses the FFTW (Fastest
	Fourier Transform in the West) external library which must be linked to MIMICA at compile time. If FFTW
	is centrally built and available on your local server, it can be directly linked to MIMICA V5 with only small
        modifications. In order to use the central FFTW library, MIMICA should be compiled with 
	Makefile_centralfft which must first be copied into your local work directory, and renamed Makefile. If the 
	environment variables FFTW_INC and FFTW_LIB are not already defined on your system, you should define those 
	variables yourself.
	
		If FFTW is not available or if you do not wish to link to a centrally
	distributed version of fftw, a local copy is provided in src/fftw-3.3.8. The local FFTW copy must be configured and 
	compiled separately. To build the local FFTW library, first cd to the src/fftw-3.3.8 directory. Then create your own 
	local build by typing successively (the whole procedure can take up to 5-10min): 

		./configure --prefix=$PWD/build --enable-mpi --enable-shared

		make 

		make install

	If FFTW had previously been compiled in this directory, you will first need to clear all existing files by typing:

		make distclear
		
		You can check that the FFTW library has been properly installed by verifying that the files 
	libfftw3.a and fftw3.f03 are present in fftw-3.3.8/build/lib and fftw-3.3.8/build/include respectively. If anything went 
	wrong during the installation, this may be due to inconsistencies in the C compiler or MPI wrapper used. Once FFTW has
	been successfully built, MIMICA can be compiled by using the default Makefile present in the main mimicav5 directory.
	Makefile can then be copied in your local work directory without any further modification. 

* 1. (b) Compilation:

		Before going any further and proceeding with MIMICA's compilation, a new environment variable must
	be defined to specify the main path to your MIMICA directory. To do this, you must simply type:

		export MIMICA=......

	where you must appropriately complete your own MIMICA path name. For future uses, you can directly add this
	line to your ./profile file.

		Compiling MIMICA is performed using a 'start' script which must be appropriately edited to compile
	with the desired physical and numerical packages. 'start' templates are provided for a variety
	of idealized test cases in the ./templates folder. The most important flags which must be edited
	in 'start' are (please refer to README.start for an exhaustive list of parameters):
                - NETCDF: path to your netcdf directory
		- NP: number of processors
		- CMPLER: choice of compiler. The default compiler is gfortran. If setenv CMPLER INTEL is uncommented,
			the intel fortran compiler is used.
		- DEBUG: if TRUE, MIMICA is compiled with debug options
		- SPMD: if TRUE, MIMICA is compiled for parallel runs (NP must then be > 1)
		- M3D: if TRUE, a 3D run is setup
		- MAXX, MAXY, MAXZ (in 2D or 3D): number of grid points in the domain along X, Y, Z
		
		Once 'start' has been copied from a template and edited, you must provide one of the two Makefiles
	available in the parent directory as described in section 1a above. Compilation is then performed by typing:
		./start
		
		All objects and modules are stored in ./build. The file compile.log created in your work directory
	contains a summary of the compilation with all the different start options selected. 
	Finally, the executable 'mimicav5_XXXXXX.exe' is created in your local MIMICA directory. XXXXXX is a flag
        corresponding to the current commit number from which your executable has been built (this enables you to keep
        track of the actual versions of your executables). 
	
		Compiling MIMICA requires various external libraries which are usually available directly on your local
	workstation/HPC system (for example with the module load command). The required libraries include a recent
	version of the compiler (intel or gcc), a compatible MPI library, Netcdf and FFTW (see above for further details). 
	The important thing to remember is to always use compatible libraries: if gfortran (gcc) is selected to compile 
	MIMICA, consistent versions of the MPI, NETCDF and FFTW libraries must be loaded.

* 2. Run:

		Running MIMICA can be done in serial by typing:

		./mimicav5_XXXXXXX.exe
		
		or in parallel with:

		mpirun -np Y mimicav5_XXXXXXX.exe	(with Y the number of processors)
		
		Please consult your HPC help page to get more information on how to properly run a parallel job.
		
		A number of input files required by MIMICA must also exist (3 in your local directory, 'start', 
        'cm.nml' and 'out.nml', as well as appropriate initial conditions in ./INCLUDE).

                The most important one is 'cm.nml' which contains a list of model options and parameters. Again,
	'cm.nml' templates are available for idealized cases in ./templates. The file must be present (copied from
        a template) in your local working directory and edited as desired. Your working directory must also contain
        an 'out.nml' file where you can define the variables you want to ouput using various keywords. An example
        is availabel in the main MIMICA repository. 
	
		Other than these 3 essential files which must be present in your working directory, you will also need 
	initial conditions for your simulation. These will be prescribed using either a specific include file, .h, which 
	must be present in ./INCLUDE, or by creating your own initial sounding, again in ./INCLUDE (the name of the initial
        sounding file is prescribed in 'cm.nml'). 

* 3. After the run:

		Output files are created and stored in ./OUTPUT. At the moment, these consist only in raw data files which
	must be further processed to produce files readable by visualization softwares. A fortran program for converting 
	the raw data into a more structured format is available in ./tools (post_process.f90). This program can also be
	used to process the 1D (horizontally averaged) profiles created by MIMICA. Another good option to process these
	raw data would be to use matlab.

		In case you need to restart MIMICA to continue the simulation, a 'restart.dat' file is created in ./INCLUDE.
	To run MIMICA starting from the state stored in 'restart.dat', there is no need to recompile MIMICA. The only thing
	required is to set the 'new_run' flag in 'cm.nml' to true (and make sure the simulation endtime is properly set). 
	
### Who do I talk to? ###

* Main developer: Julien Savre (LMU Munich - julien.savre@lmu.de)
* Code administrator: Matthias Brakebusch (ACES Stockholm University - matthias.brakebusch@aces.su.se)
* Other main contributors: Ben Murphy (ben.n.murphy@gmail.com), Annica Ekman (annica@misu.su.se)
