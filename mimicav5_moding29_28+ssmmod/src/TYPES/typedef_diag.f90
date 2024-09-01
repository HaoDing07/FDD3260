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
!  TYPEDEF_DIAG.F:                   
!
!  Purpose:
!	Typedef for diagnostics	: The diagnostic quantities stored in
!	in this structure correspond to the tendencies used to advance
!	the prognostic quantities as well as possible derived varibles.		  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!		
! ================================================================
	
	module typedef_diag
	
	USE gridno
	USE shared_data
	
	type :: gradient
		! ---
		real :: u	 ! u wind
		real :: v	 ! v wind
		real :: w	 ! w wind
		real :: qt	 ! Total water mixing ratio
		real :: pt	 ! ice/liquid potential temperature
	end type gradient
	
	type :: microdiag
		! ---
		real, dimension(:,:,:), allocatable :: q
		real, dimension(:,:,:), allocatable :: n
	end type microdiag
	
	type :: aerodiag
		! ---
		real, dimension(:,:,:), allocatable :: n
		real, dimension(:,:,:), allocatable :: m
	end type aerodiag
	
	type :: tendencies
		! ---
		real, dimension(:,:,:), allocatable :: u       ! u wind
		real, dimension(:,:,:), allocatable :: v       ! v wind
		real, dimension(:,:,:), allocatable :: w       ! w wind
		real, dimension(:,:,:), allocatable :: qt      ! Total water mixing ratio
		real, dimension(:,:,:), allocatable :: pt      ! potential temperature
		real, dimension(:,:,:), allocatable :: ptv     ! Virtual potential temperature
		real, dimension(:,:,:), allocatable :: p       ! Pressure
		real, dimension(:,:,:), allocatable :: k       ! Kinetic energy
		real, dimension(:,:,:), allocatable :: qc      ! Droplet mixing ratio
		real, dimension(:,:,:), allocatable :: qr      ! Rain mixing ratio
		real, dimension(:,:,:), allocatable :: qi      ! Ice mixing ratio
		real, dimension(:,:,:,:), allocatable :: sca   ! Other scalars
		type(gradient), dimension(:,:,:), allocatable :: grad
		type(microdiag), dimension(:), allocatable :: micro
		type(aerodiag), dimension(:), allocatable :: aero
	end type tendencies
	
	type :: diagnosvar
		! ---
		real, dimension(:,:,:), allocatable :: z	! Radar reflectivity
		real, dimension(:,:,:), allocatable :: pb	! Buoyant pressure pert.
		real, dimension(:,:,:), allocatable :: pd	! Dynamical pressure pert.
		real, dimension(:,:,:), allocatable :: pnh	! Non hydrostatic pressure pert.
		real, dimension(:,:,:), allocatable :: beff	! Effective buoyancy
		real, dimension(:,:,:), allocatable :: div, qconv ! Horizontal divergence
		real, dimension(:,:,:), allocatable :: vortx, vorty, vortz	! Vorticity
                real, dimension(:,:,:), allocatable :: mf, cc, ccm, mfm   ! Convective cloud diag.
        end type diagnosvar

        type :: entrainment
		real, dimension(:,:,:), allocatable :: dent 
		real, dimension(:,:,:), allocatable :: ddet 
		real, dimension(:,:,:), allocatable :: activ    ! Activity
		real, dimension(:,:,:), allocatable :: went	! Entrainment rate
		real, dimension(:,:,:), allocatable :: went1	! Entrainment rate
		real, dimension(:,:,:), allocatable :: went2	! Entrainment rate
		real, dimension(:,:,:), allocatable :: went3	! Entrainment rate
		real, dimension(:,:,:), allocatable :: went4	! Entrainment rate
		real, dimension(:,:,:), allocatable :: went5	! Entrainment rate
		real, dimension(:,:,:), allocatable :: cent	! Entrainment rate
		real, dimension(:,:,:), allocatable :: wdet	! Detrainment rate
		real, dimension(:,:,:), allocatable :: wdet1	! Detrainment rate
		real, dimension(:,:,:), allocatable :: wdet2	! Detrainment rate
		real, dimension(:,:,:), allocatable :: wdet3	! Detrainment rate
		real, dimension(:,:,:), allocatable :: wdet4	! Detrainment rate
		real, dimension(:,:,:), allocatable :: wdet5	! Detrainment rate
		real, dimension(:,:,:), allocatable :: cdet	! Detrainment rate
        end type entrainment

	type :: filtered
	  	real :: epse
	        real :: qte
		real :: pte
		real :: ptve
		real :: buoye
		real :: beffe
		real :: we
		real :: epsd
	        real :: qtd
		real :: ptd
		real :: ptvd
		real :: buoyd
		real :: beffd
		real :: wd
		real :: epsa
	        real :: qta
		real :: pta
		real :: ptva
		real :: buoya
		real :: beffa
		real :: wa
	end type filtered
	
	type :: plumes
	        real :: qt, ptv, w
	        real :: qte, ptve, we
	        real :: qtb, ptvb, wb
		real :: cov
		real :: mfl
		real :: dadt
		real :: dmfl
		real :: mu
		real :: eq
		real :: et
		real :: ew
		real :: eqs
		real :: ets
		real :: ews	
		real :: alphaq
		real :: alphat
		real :: alphaw
		real :: alphab
		real :: alphap
		real :: beta
		real :: gamma
		real :: delta1, delta2
		real :: wvar, wske
		real :: ri, mri
	end type plumes
	
	interface alloc
		module procedure allocate_diag, allocate_diagnos, allocate_entrain
	end interface
	
	interface dealloc
		module procedure deallocate_diag, deallocate_diagnos, deallocate_entrain
	end interface

	CONTAINS
	
	subroutine allocate_diag ( diag )
	      type(tendencies), dimension(ndiag), intent(inout) :: diag
	      integer :: k, l
	
	      do k = 1, ndiag
	        if (out_diagu) allocate ( diag(k)%u(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_diagv) allocate ( diag(k)%v(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_diagw) allocate ( diag(k)%w(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_diagq) allocate ( diag(k)%qt(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_diagt) allocate ( diag(k)%pt(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_diagtv) allocate ( diag(k)%ptv(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_diagp) allocate ( diag(k)%p(ip_start:ip_end,jp_start:jp_end,1:nz) ) 
	        if (out_diagk) allocate ( diag(k)%k(ip_start:ip_end,jp_start:jp_end,1:nz) ) 
	        if (out_diagl.and.lmicro > 0) allocate ( diag(k)%qc(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_diagr.and.lmicro > 0) allocate ( diag(k)%qr(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_diagi.and.lmicro > 1) allocate ( diag(k)%qi(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_diags.and.nscal > 0) allocate ( diag(k)%sca(ip_start:ip_end,jp_start:jp_end,1:nz,1:nscal) )
	        if (out_grad) allocate ( diag(k)%grad(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        
	        if (out_micro.and.lmicro>0) then
	          allocate ( diag(k)%micro(1:nhydro) )
	          do l = 1, nhydro
	            allocate ( diag(k)%micro(l)%q(ip_start:ip_end,jp_start:jp_end,1:nz) )
	            allocate ( diag(k)%micro(l)%n(ip_start:ip_end,jp_start:jp_end,1:nz) )
	          enddo
	        endif
	        
	        if (out_aero.and.out_diaga) then
	          allocate ( diag(k)%aero(1:nmode) )
	          do l = 1, nmode
	            allocate ( diag(k)%aero(l)%n(ip_start:ip_end,jp_start:jp_end,1:nz) )
	            allocate ( diag(k)%aero(l)%m(ip_start:ip_end,jp_start:jp_end,1:nz) )
	          enddo
	        endif
	      enddo

	return
	end subroutine allocate_diag
	
	subroutine deallocate_diag ( diag )
	      type(tendencies), dimension(ndiag), intent(inout) :: diag
	      integer :: k, l
	
	      do k = 1, ndiag
	        if (out_diagu) deallocate ( diag(k)%u )
	        if (out_diagv) deallocate ( diag(k)%v )
	        if (out_diagw) deallocate ( diag(k)%w )
	        if (out_diagq) deallocate ( diag(k)%qt )
	        if (out_diagt) deallocate ( diag(k)%pt )
	        if (out_diagtv) deallocate ( diag(k)%ptv )
	        if (out_diagp) deallocate ( diag(k)%p )
	        if (out_diagk) deallocate ( diag(k)%k )
	        if (out_diagl.and.lmicro > 0) deallocate ( diag(k)%qc )
	        if (out_diagr.and.lmicro > 0) deallocate ( diag(k)%qr )
	        if (out_diagi.and.lmicro > 1) deallocate ( diag(k)%qi )
	        if (out_diags.and.nscal > 0) deallocate ( diag(k)%sca )
	        if (out_grad) deallocate ( diag(k)%grad )
	        
	        if (out_micro.and.lmicro>0) then
	          do l = 1, nhydro
	            deallocate ( diag(k)%micro(l)%q )
	            deallocate ( diag(k)%micro(l)%n )
		  enddo
	          deallocate ( diag(k)%micro )
	        endif
	        
	        if (out_aero.and.out_diaga) then
	          do l = 1, nmode
	            deallocate ( diag(k)%aero(l)%n )
	            deallocate ( diag(k)%aero(l)%m )
		  enddo
	          deallocate ( diag(k)%aero )
	        endif
	      enddo

	return
	end subroutine deallocate_diag
	
	subroutine allocate_diagnos ( diag )
	      type(diagnosvar), intent(inout) :: diag
	
	      if ( spec_diag ) then
	        if (out_beff) allocate ( diag%beff(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_dp) allocate ( diag%pb(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_dp) allocate ( diag%pd(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_dp) allocate ( diag%pnh(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_z) allocate ( diag%z(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_div) allocate ( diag%div(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_div) allocate ( diag%qconv(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_vort) allocate ( diag%vortx(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_vort) allocate ( diag%vorty(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_vort) allocate ( diag%vortz(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_mf) allocate ( diag%mf(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_mf) allocate ( diag%cc(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_mf) allocate ( diag%mfm(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_mf) allocate ( diag%ccm(ip_start:ip_end,jp_start:jp_end,1:nz) )
!
	        if (out_beff) diag%beff = 0. 
	        if (out_dp) diag%pb = 0. 
	        if (out_dp) diag%pd = 0.
	        if (out_dp) diag%pnh = 0.
	        if (out_z) diag%z = 0. 
	        if (out_div) diag%div = 0.
	        if (out_div) diag%qconv = 0.
	        if (out_vort) diag%vortx = 0.
	        if (out_vort) diag%vorty = 0.
	        if (out_vort) diag%vortz = 0.
	        if (out_mf) diag%mf = 0.
	        if (out_mf) diag%cc = 0.
	        if (out_mf) diag%mfm = 0.
	        if (out_mf) diag%ccm = 0.
	      endif
	      
	return
	end subroutine allocate_diagnos
	
        subroutine allocate_entrain ( diag )
              type(entrainment), intent(inout) :: diag

              if ( spec_diag ) then
	        if (out_ent) allocate ( diag%activ(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%dent(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%ddet(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%cent(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%cdet(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%went(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%wdet(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%went1(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%went2(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%went3(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%went4(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%went5(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%wdet1(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%wdet2(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%wdet3(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%wdet4(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_ent) allocate ( diag%wdet5(ip_start:ip_end,jp_start:jp_end,1:nz) )
!
                if (out_ent) diag%activ = 0.
	        if (out_ent) diag%dent = 0. 
	        if (out_ent) diag%ddet = 0.
	        if (out_ent) diag%went = 0.
	        if (out_ent) diag%went1 = 0.
	        if (out_ent) diag%went2 = 0. 
	        if (out_ent) diag%went3 = 0.
	        if (out_ent) diag%went4 = 0. 
	        if (out_ent) diag%went5 = 0. 
	        if (out_ent) diag%cent = 0. 
	        if (out_ent) diag%wdet = 0.
	        if (out_ent) diag%wdet1 = 0.
	        if (out_ent) diag%wdet2 = 0. 
	        if (out_ent) diag%wdet3 = 0.
	        if (out_ent) diag%wdet4 = 0. 
	        if (out_ent) diag%wdet5 = 0. 
	        if (out_ent) diag%cdet = 0.
              endif

        return
        end subroutine allocate_entrain

	subroutine deallocate_diagnos ( diag )
	      type(diagnosvar), intent(inout) :: diag
	
	      if ( spec_diag ) then
	        if (out_beff) deallocate ( diag%beff )
	        if (out_dp) deallocate ( diag%pb )
	        if (out_dp) deallocate ( diag%pd )
	        if (out_dp) deallocate ( diag%pnh )
	        if (out_z) deallocate ( diag%z )
	        if (out_div) deallocate ( diag%div )
	        if (out_div) deallocate ( diag%qconv )
	        if (out_vort) deallocate ( diag%vortx )
	        if (out_vort) deallocate ( diag%vorty )
	        if (out_vort) deallocate ( diag%vortz )
	        if (out_mf) deallocate ( diag%mf )
	        if (out_mf) deallocate ( diag%cc )
	        if (out_mf) deallocate ( diag%mfm )
	        if (out_mf) deallocate ( diag%ccm )
	      endif
	      
	return
	end subroutine deallocate_diagnos

        subroutine deallocate_entrain ( diag )
	      type(entrainment), intent(inout) :: diag

              if ( spec_diag ) then
	        if (out_ent) deallocate ( diag%activ )
                if (out_ent) deallocate ( diag%dent )
	        if (out_ent) deallocate ( diag%ddet )
	        if (out_ent) deallocate ( diag%went )
	        if (out_ent) deallocate ( diag%went1 )
	        if (out_ent) deallocate ( diag%went2 )
	        if (out_ent) deallocate ( diag%went3 )
	        if (out_ent) deallocate ( diag%went4 )
	        if (out_ent) deallocate ( diag%went5 )
	        if (out_ent) deallocate ( diag%cent )
	        if (out_ent) deallocate ( diag%wdet )
	        if (out_ent) deallocate ( diag%wdet1 )
	        if (out_ent) deallocate ( diag%wdet2 )
	        if (out_ent) deallocate ( diag%wdet3 )
	        if (out_ent) deallocate ( diag%wdet4 )
	        if (out_ent) deallocate ( diag%wdet5 )
	        if (out_ent) deallocate ( diag%cdet )
              endif

        return
        end subroutine deallocate_entrain

	end module typedef_diag

