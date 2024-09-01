!
!------------------------------------------------------------------------!
! This file is part of MIMICA                                            !
!                                                                        !
! Copyright 2017-2021 Julien Savre                                       ! 
!                                                                        !
! This program is free software: you can redistribute it and/or modify   !
! it under the terms of the GNU General Public License as published by   !
! the Free Software Foundation, either version 3 of the License, or      !
! (at your option) any later version.                                    !
!                                                                        !
! This program is distributed in the hope that it will be useful,        !
! but WITHOUT ANY WARRANTY; without even the implied warranty of         !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          !
! GNU General Public License for more details.                           !
!                                                                        !
! You should have received a copy of the GNU General Public License      !
! along with this program.  If not, see <http://www.gnu.org/licenses/>.  !
!------------------------------------------------------------------------!
!

#include "ctrparam.h"

! ================================================================
!
!  TYPEDEF_AEROSOL.F:                   
!
!  Purpose:
!	Typedef for derived-type data: aerosol			  
!
!  Author
!	Annica Ekman based on module by Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!		
! ================================================================
	
  module typedef_aerosol_new

  use gridno
  use shared_data
!
!  Dimensions
!
	integer, parameter :: nelem = 4          !> number of aerosol species used in MIMICA
	integer, parameter :: nmo = 25, npress = 54, nrh = 46
		
	type :: elementary
	  real :: kappa               !> hygroscopicity
	  real :: rho                 !> density
	  real :: mw                  !> molecular weight
	  character(len=10) :: name
	end type elementary

	type :: aero_init0
		! Type to define inital properties of a aerosol mode 
	    logical, dimension(nelem) :: present         !> bool list; defined is species is present in mode
	    real, dimension(nelem)    :: frac            !> defines fraction of a species in mode
		real    		      :: n0                  !> initial number concentration of mode
		real                  :: n_sfc_source        !> suface aerosol number source (from cm.nml) 
	end type aero_init0

	type :: aero_size
		! Type to define properties of lognormal PDF
		real :: rmean            !> Initial geometric mean radius 
		real :: mmean            !> Initial mass concentration
		real :: sigma            !> Geometric standard deviation
	end type aero_size
	
	type :: aero_init
		! Type to define aerosol mode
	    integer	      :: nelem
	    type(aero_init0)  :: init
	    type(aero_size)   :: size
	end type aero_init
	
	type :: aeroprop
	    ! Type for properties of aerosol mode
	    integer	      :: nelem
	    type(aero_init0)  :: init
	    type(aero_size)   :: size
	    type(elementary)  :: mix
	    real 	      :: depv0
	end type aeroprop

	type :: aero_1d
	    real, dimension(nz) :: n, m, ma       !> number, mass, and activated mass of aerosol
	end type aero_1d

	type :: aero_3d
	    real, allocatable, dimension(:,:,:) :: n, m, ma  !> number, mass, and activated mass of aerosol
	end type aero_3d 

	type :: aero_param
		real :: critn
		real :: nucint
		real :: aitint
		real :: bcagemass
	end type aero_param

	type :: aero_bulk
		real :: mmol
		real :: rho
		real :: nuphi
		real :: beta
		real :: b
	end type aero_bulk

	!
	! Flags for aerosol processes
	!
	type aero_flg_list
		logical :: tran        !< transport of aerosol on/off
		logical :: act         !< activation of aerosol on/off
		logical :: act_scv     !< activation scavenging of aerosol on/off
		logical :: imp_scv     !< impaction scavenging of aerosol on/off
		logical :: reg         !< regeneration of aerosol on/off
		logical :: chem        !< aerosol chemisty on/off
		logical :: any         !< flag to sumarize if any aerosol process is sumulated
	end type aero_flg_list

	
	interface assignment (=)
		module procedure aero_to_aero_init, aero_to_aero_size, aero_to_aero_3d
	end interface
		
        interface operator (+)
        module procedure aero_plus_aero_3d
        end interface
	
	interface alloc
		module procedure allocate_aero
	end interface
	
	interface dealloc
		module procedure deallocate_aero
	end interface
	
	interface reset
		module procedure reset_aero3d
	end interface

	CONTAINS
	
	subroutine aero_to_aero_init (aero1, aero2)
	      type (aero_init0), intent(out)	 :: aero1
	      type (aero_init0), intent(in)	 :: aero2
	
	      aero1%present(:)    = aero2%present(:)
	      aero1%frac(:)       = aero2%frac(:)
	      aero1%n0            = aero2%n0
	      aero1%n_sfc_source  = aero2%n_sfc_source
	return
	end subroutine aero_to_aero_init
	
	subroutine aero_to_aero_size (aero1, aero2)
	      type (aero_size), intent(out)	 :: aero1
	      type (aero_size), intent(in)	 :: aero2
	
	      aero1%rmean         = aero2%rmean
	      aero1%mmean         = aero2%mmean
	      aero1%sigma         = aero2%sigma
	return
	end subroutine aero_to_aero_size
	
	subroutine aero_to_aero_3d (aero1, aero2)
	      integer :: im, iel
	      type (aero_3d), intent(out) :: aero1(nmode)
	      type (aero_3d), intent(in)  :: aero2(nmode)
	
	      do im = 1, nmode            
	        aero1(im)%n = aero2(im)%n
	        aero1(im)%m = aero2(im)%m
	        aero1(im)%ma = aero2(im)%ma
              enddo
	return
	end subroutine aero_to_aero_3d

        function  aero_plus_aero_3d (aero1, aero2) result(aero3)
              integer :: im, iel
              type (aero_3d), intent(in) :: aero1(nmode)
              type (aero_3d), intent(in)  :: aero2(nmode)
              type (aero_3d)  :: aero3(nmode)

              do im = 1, nmode
                aero3(im)%n = aero1(im)%n + aero2(im)%n
                aero3(im)%m = aero1(im)%m + aero2(im)%m
                aero3(im)%ma = aero1(im)%ma + aero2(im)%ma
              enddo
        return
        end function aero_plus_aero_3d

! ============================	
  subroutine add_to_aero3d (aerotmp, q, c1, im, iel)
! ============================	
  
  integer :: i, j, k, im
  integer, optional :: iel
  character(len=1) :: c1
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: q
  type (aero_3d), dimension(nmode) :: aerotmp     

  do k = 1, nz
    do j = 1, jp_start, jp_end
      do i = ip_start, ip_end
        if (c1=='n') aerotmp(im)%n(i,j,k) = aerotmp(im)%n(i,j,k) + q(i,j,k)
        if (c1=='m') aerotmp(im)%m(i,j,k) = aerotmp(im)%m(i,j,k) + q(i,j,k)
        if (c1=='a') aerotmp(im)%ma(i,j,k) = aerotmp(im)%ma(i,j,k) + q(i,j,k)
      enddo
    enddo
  enddo

  end subroutine add_to_aero3d   

! ============================	
  subroutine allocate_aero (aerotmp)
! ============================	

  	type(aero_3d), dimension(1:nmode), intent(inout) :: aerotmp
	integer :: h

	do h = 1, nmode
	  allocate ( aerotmp(h)%n(ip_start:ip_end,jp_start:jp_end,1:nz) )
	  allocate ( aerotmp(h)%m(ip_start:ip_end,jp_start:jp_end,1:nz) )
	  allocate ( aerotmp(h)%ma(ip_start:ip_end,jp_start:jp_end,1:nz) )

	  aerotmp(h)%n = 0.
	  aerotmp(h)%m = 0.
	  aerotmp(h)%ma = 0.
        enddo

  return
  end subroutine allocate_aero

! ============================	
  subroutine deallocate_aero (aerotmp)
! ============================	

  	type(aero_3d), dimension(1:nmode), intent(inout) :: aerotmp
	integer :: h

	do h = 1, nmode
	  deallocate ( aerotmp(h)%n )
	  deallocate ( aerotmp(h)%m )
	  deallocate ( aerotmp(h)%ma )
        enddo

  return
  end subroutine deallocate_aero

! ============================	
  subroutine reset_aero3d (aerotmp)
! ============================	
  
  integer :: i, j, k, im
  type (aero_3d), dimension(nmode) :: aerotmp   
  
  do k = 1, nz
    do j = jp_start, jp_end
      do i = ip_start, ip_end
        do im = 1, nmode
          aerotmp(im)%n(i,j,k) = 0.
          aerotmp(im)%m(i,j,k) = 0.
          aerotmp(im)%ma(i,j,k) = 0.
        enddo  
      enddo
    enddo
  enddo
  
  return
  end subroutine
  	  	    			
  end module typedef_aerosol_new

