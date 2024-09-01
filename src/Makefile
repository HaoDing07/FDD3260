AERODIR=./AEROCHEM
INTEGDIR=./INTEG
IODIR=./IO
MAINDIR=./MAIN
MICRODIR=./MICRO
PHYSDIR=./PHYS
PRESSUREDIR=./PRESSURE
FFTDIR=./fftmpi/src
TOOLSDIR=./TOOLS
TYPESDIR=./TYPES

all: create-fft create-types create-tools create-aero create-micro create-phys create-pressure\
	create-io create-integ create-main

create-aero:
	cd $(AERODIR) ; $(MAKE)

create-integ:
	cd $(INTEGDIR) ; $(MAKE)

create-io:
	cd $(IODIR) ; $(MAKE)

create-main:
	cd $(MAINDIR) ; $(MAKE)

create-micro:
	cd $(MICRODIR) ; $(MAKE)

create-phys:
	cd $(PHYSDIR) ; $(MAKE)

create-tools:
	cd $(TOOLSDIR) ; $(MAKE)

create-fft:
	cd $(FFTDIR) ; $(MAKE)

create-pressure:
	cd $(PRESSUREDIR) ; $(MAKE)

create-types:
	cd $(TYPESDIR) ; $(MAKE)

clean: clean-fft clean-types clean-aero clean-integ clean-io clean-main clean-micro\
	clean-phys clean-pressure clean-tools
	
clean-aero:
	cd $(AERODIR) ; $(MAKE) clean-aero
	
clean-integ:
	cd $(INTEGDIR) ; $(MAKE) clean-integ

clean-io:
	cd $(IODIR) ; $(MAKE) clean-io

clean-main:
	cd $(MAINDIR) ; $(MAKE) clean-main

clean-micro:
	cd $(MICRODIR) ; $(MAKE) clean-micro

clean-phys:
	cd $(PHYSDIR) ; $(MAKE) clean-phys

clean-fft:
	cd $(FFTDIR) ; $(MAKE) clean-fft

clean-pressure:
	cd $(PRESSUREDIR) ; $(MAKE) clean-pressure

clean-tools:
	cd $(TOOLSDIR) ; $(MAKE) clean-tools

clean-types:
	cd $(TYPESDIR) ; $(MAKE) clean-types
