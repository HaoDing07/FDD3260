#!/bin/csh -f
echo "Compile gridno.f90"

cat >! $MIMICA/src/TYPES/gridno.f90 << EOF
#include "ctrparam.h"

      module gridno
      
! --- from GRIDNO.H:
! ---  A header file contains parameter statements 
! ---             for Grid Numbers
! ---       cpp F90/F95 version

      integer, parameter :: nprocx = $NPX
      integer, parameter :: nprocy = $NPY

      integer, parameter :: nx = $MAXX+5
#ifdef MODEL_3D
      integer, parameter :: ny = $MAXY+5
#else
      integer, parameter :: ny = 1
#endif
      integer, parameter :: nz = $MAXZ, nzi =nz-1
 
      character(len = 100) :: datadir = '$DATADIR'
      character(len = 100) :: incdir = '$INCDIR'
      character(len = 10 ) :: version = '$VERSION'
      
      integer :: ip_start         ! begining i of a patch including ghost cells
      integer :: ip_end           ! end i      of a patch including ghost cells
      integer :: jp_start         ! begining j of a patch including ghost cells
      integer :: jp_end           ! end j      of a patch including ghost cells

      integer :: it_start         ! begining i of a tile
      integer :: it_end           ! end i      of a tile
      integer :: jt_start         ! begining j of a tile
      integer :: jt_end           ! end j      of a tile
                  
      integer  :: nx_p
      integer  :: nx_t
      integer  :: ny_t
      integer  :: ny_p
       
      integer, parameter :: nscal = $NSCAL	! Number of passive scalars
      integer, parameter :: npart = $NPART	! Number of passive scalars

      integer, parameter :: nhab = 9		! Number of ice habits
      integer, parameter :: ndiag  = 10 	! Number of diagnostic variables
      integer, parameter :: ntsout = 32		! Number of time-series variables
      integer, parameter :: ntmpout= 7		! Number of user defined output variables
      integer, parameter :: nmark = 0		! Number of special markers for lagrangian parcels

      integer, parameter :: n_trop = 500

      integer :: drop            ! Index of cloud drops 
        parameter(drop=1)
      integer :: rain            ! Index of rain drops
        parameter(rain=2)      
      integer :: ice             ! Index of ice particles
        parameter(ice=3)      
      integer :: grau            ! Index of graupel particles
        parameter(grau=4)      
      integer :: snow            ! Index of snow particles
        parameter(snow=5)      
      integer :: hail            ! Index of hail particles
        parameter(hail=6)      

      end module gridno
EOF
