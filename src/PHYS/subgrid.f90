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

!  =================================================
!
!	SUBGRID.F:  A package fturbutor calculating 
!		    eddy coefficient - K and subgrid
!		    scale contributions of winds as well
!		    as scalars, Using 4th order differences
!
!	Author:   Julien Savre, Misu
! 
!  ==================================================

  module subgrid
  
  USE shared_all
  USE advection
  USE surfacemod
  USE allocation
  USE boundary_conditions
  USE gradients
  USE sources
  USE allocation
!
  IMPLICIT NONE
  
  private

  public :: subgrid_interface, calc_smago, tkeinit, efud, eff0, effs

  CONTAINS

! ===============================================
  subroutine subgrid_interface ( pressurel, windl, statel, hydromtrl )
! ===============================================
!
  type (atm_pressure) :: pressurel
  type (atm_winds) :: windl
  type (atm_state) :: statel
  type (hydrometeor), dimension(1:nhydro) :: hydromtrl
!
  integer :: h
!
  if (verbose > 0) call write_debug('Starting subgrid_interface')
!
!----------------------------------------------------------!
!                Calculate eddy viscosity                  !
!----------------------------------------------------------!
!
  call calc_subgrid ( pressurel, windl, statel, hydromtrl )
!
!-----------------------------------------------------------!
!              	        Surface fluxes                      !
!-----------------------------------------------------------!
!
  if ( isurf >= 0 ) call calc_surface ( pressurel, windl, statel, hydromtrl )
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating subgrid_interface')
!
return
end
!
! ===============================================
  subroutine calc_subgrid ( pressurel, windl, statel, hydromtrl )
! ===============================================
!
  type(atm_pressure), intent(in) :: pressurel
  type(atm_winds), intent(in) :: windl
  type(atm_state), intent(in) :: statel
  type(hydrometeor), dimension(:), intent(in) :: hydromtrl
!
  real, dimension(:,:,:), allocatable :: N2, S2
  real, dimension(:,:,:,:,:), allocatable :: dij
!
  call alloc ( dij, 3, 3 )
  call alloc ( N2 )
  call alloc ( S2 )
!
  if (verbose > 1) call write_debug('Starting calc_subgrid')
!
!----------------------------------------------------------!
!	 Calculate and store stability parameters          !
!----------------------------------------------------------!
!
  if (diff == 0.) then
!
    call cal_dij ( windl, dij=dij ) 
!
    call stability ( pressurel, statel, hydromtrl, dij, N2=N2, S2=S2 )
!
!  Calculate eddy viscosity 
!
#ifdef TKE
    call tkecomp ( pressurel, windl, statel, hydromtrl, N2, S2 ) 
!
    call calc_tke ( N2 )
#else
    call calc_smago ( N2, S2 )
#endif
!
!  Diagnostics
!
    if ( out_tke ) then
      turbu_diag%N2 = N2
      turbu_diag%S2 = S2
    endif
!
!  Fixed diffusion coefficient
!
  else
    turbu%fkv = abs(diff)
  endif
!
  call dealloc ( dij )
  call dealloc ( N2 )
  call dealloc ( S2 )
!
  if (verbose > 1) call write_debug('Terminate calc_subgrid')
!
return
end
!
! ==================================================
  subroutine tkeinit ( pressurel, windl, statel, hydromtrl )
! ==================================================
!
  type (atm_pressure), intent(in) :: pressurel
  type (atm_winds), intent(in) :: windl
  type (atm_state), intent(in) :: statel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
!
  integer :: k
  real, dimension(:,:,:,:,:), allocatable :: dij
  real, dimension(:,:,:), allocatable :: N2, S2
!
!----------------------------------------------------------!
!          	  Calculate Initial ksgs                   !
!----------------------------------------------------------!
!
!  Allocate
!
  call alloc ( dij, 3, 3 )
  call alloc ( N2 )
  call alloc ( S2 )
!
  call cal_dij ( windl, dij=dij )
!
  call stability ( pressurel, statel, hydromtrl, dij, N2=N2, S2=S2 )
!
  call calc_smago ( N2, S2 )
!
!  Deallocate
!
  call dealloc ( dij )
  call dealloc ( N2 )
  call dealloc ( S2 )
!
!-----------------------------------------------------------!
!
return
end
!
!  =======================================
   subroutine efud ( pressurel, windl, statel, hydromtrl )
!  =======================================
!
! ------------------------------------------------------
! --- Subroutine for calculating eddy terms of vectors
! ------------------------------------------------------
!
  type (atm_pressure), intent(in) :: pressurel
  type (atm_winds), intent(in) :: windl
  type (atm_state), intent(in) :: statel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
!
  integer :: i, j, k, i1, i2, j1, j2
  real  :: zz, dx1, dy1, dt1, func, f_an, windm
!
  real, allocatable, dimension(:,:,:) :: ev, eh, flv, fkv, fkh
  real, allocatable, dimension(:,:,:,:,:) :: dij, Flux
!
  if (verbose > 1) call write_debug('Enter efud')
!
!----------------------------------------------------------!
!                 Calculate surface fluxes                 !
!----------------------------------------------------------!
!
!  Allocate
!
  call alloc ( dij, 3, 3 )
  call alloc ( Flux, 3, 3 )
  call alloc ( eh )
  call alloc ( ev )
  call alloc ( flv )
  call alloc ( fkv )
  call alloc ( fkh )
!
!  Calculate stress tensor 
!
  call cal_dij ( windl, dij=dij ) 
!
!  Non-linear closure
!
  if (nl_sgs) then
    call non_linear ( pressurel, windl, statel, hydromtrl, dij, fkv ) 
  else
    fkv = pressurel%dens*turbu%fkv
    dij = 2.*dij
  endif
!
!  Eddy viscosity damping
!
  if ( pran < 0.0 ) then
    do k = 1, nz
      func = max(1. - exp(6.*(z0(k) - zdec)/zdec),0.)
      fkv(:,:,k) = func*fkv(:,:,k)
    enddo
  endif
!
!  anisotropy factor  
!
  fkh = fkv
  if (anis_k) then
    do k = 1, nz
      f_an = ( (dx*dy)**(1./2.) / (dz*fdz0(k)) )**2.
      fkh(:,:,k) = f_an*fkh(:,:,k)
    enddo
  endif
!
!----------------------------------------------------------!
!     	       Calculating diffusive fluxes    	     	   !
!----------------------------------------------------------!
!
!  Indices
!
  i1=it_start-1
  i2=it_end
  j1=1
  j2=1
#ifdef MODEL_3D
  j1=jt_start-1
  j2=jt_end
#endif
!
  dx1  = 1./dx
  dy1  = 1./dy  
  dt1  = 1./dt
!
!  For u velocity
!
  do k=1,nz-1
    do j=j1,j2
      do i=i1,i2
        Flux(i+1,j,k,1,1) = fkh(i,j,k)*dij(i,j,k,1,1)
!      
#ifdef MODEL_3D
        Flux(i,j+1,k,1,2) = 0.25*(fkh(i,j,k)*dij(i,j,k,1,2) + fkh(i,j+1,k)*dij(i,j+1,k,1,2) + fkh(i-1,j,k)*dij(i-1,j,k,1,2) + fkh(i-1,j+1,k)*dij(i-1,j+1,k,1,2)) 
#endif
!
        Flux(i,j,k+1,1,3) = 0.25*(fkv(i,j,k)*dij(i,j,k,1,3) + fkv(i-1,j,k)*dij(i-1,j,k,1,3) + fkv(i-1,j,k+1)*dij(i-1,j,k+1,1,3) + fkv(i,j,k+1)*dij(i,j,k+1,1,3))
      end do
    end do
  end do
!
!  For v velocity
!
#ifdef MODEL_3D
  do k=1,nz-1
    do j=j1,j2
      do i=i1,i2
        Flux(i+1,j,k,2,1) = 0.25*(fkh(i,j,k)*dij(i,j,k,2,1) + fkh(i,j-1,k)*dij(i,j-1,k,2,1) + fkh(i+1,j-1,k)*dij(i+1,j-1,k,2,1) + fkh(i+1,j,k)*dij(i+1,j,k,2,1))
!
        Flux(i,j+1,k,2,2) = fkh(i,j,k)*dij(i,j,k,2,2)
!      
        Flux(i,j,k+1,2,3) = 0.25*(fkv(i,j,k)*dij(i,j,k,2,3) + fkv(i,j-1,k)*dij(i,j-1,k,2,3) + fkv(i,j-1,k+1)*dij(i,j-1,k+1,2,3) + fkv(i,j,k+1)*dij(i,j,k+1,2,3))
      end do
    end do
  end do
#endif
!
!  For w velocity (all fluxes to 0. at the surface, k = 1)
!
  do k=2,nz-1
    do j=j1,j2
      do i=i1,i2
        Flux(i+1,j,k,3,1) = 0.25*(fkh(i,j,k)*dij(i,j,k,3,1) + fkh(i,j,k-1)*dij(i,j,k-1,3,1) + fkh(i+1,j,k-1)*dij(i+1,j,k-1,3,1) + fkh(i+1,j,k)*dij(i+1,j,k,3,1))
!      
#ifdef MODEL_3D
        Flux(i,j+1,k,3,2) = 0.25*(fkh(i,j,k)*dij(i,j,k,3,2) + fkh(i,j,k-1)*dij(i,j,k-1,3,2) + fkh(i,j+1,k-1)*dij(i,j+1,k-1,3,2) + fkh(i,j+1,k)*dij(i,j+1,k,3,2))
#endif
!
        Flux(i,j,k+1,3,3) = fkv(i,j,k)*dij(i,j,k,3,3)
      end do
    end do
  end do
  Flux(i1:i2,j1:j2,2,3,3) = fkv(i1:i2,j1:j2,1)*dij(i1:i2,j1:j2,1,3,3)
!
!  Surface Fluxes
!
  if (momsf > 0) then
    do j=j1,j2
      do i=i1,i2
        windm = (0.5*(windl%u(i+1,j,1)+windl%u(i,j,1)))**2.
#ifdef MODEL_3D
        windm = windm + (0.5*(windl%v(i,j+1,1)+windl%v(i,j,1)))**2.
#endif
        windm = max(min_w,sqrt(windm))
!
        Flux(i,j,1,1,3) = surf%mflux(i,j) * 0.5*(windl%u(i+1,j,1)+windl%u(i,j,1)) / windm
#ifdef MODEL_3D
        Flux(i,j,1,2,3) = surf%mflux(i,j) * 0.5*(windl%v(i,j+1,1)+windl%v(i,j,1)) / windm
#endif
      end do
    end do
  endif
!
#ifdef CHANNEL
  if (jp_start == 1) then
    Flux(:,jt_start,:,1,2) = 0.
    Flux(:,jt_start,:,3,2) = 0.
  endif
!
  if (jp_end == ny) then
    Flux(:,jt_end+1,:,1,2) = 0.
    Flux(:,jt_end+1,:,3,2) = 0.
  endif
#endif
!
!----------------------------------------------------------!
!      		   Assemble SGS fluxes    		   !
!----------------------------------------------------------!
!
!  U velocity
!
  eh = 0.0
  ev = 0.0
  flv = 0.0
  do k=1,nz-1
    do j=jt_start,jt_end
      do i=it_start,it_end
	eh(i,j,k) = (Flux(i+1,j,k,1,1) - Flux(i,j,k,1,1))*dx1
!
#ifdef MODEL_3D
	eh(i,j,k) = eh(i,j,k) + (Flux(i,j+1,k,1,2) - Flux(i,j,k,1,2))*dy1
#endif
!
	ev(i,j,k) = (Flux(i,j,k+1,1,3) - Flux(i,j,k,1,3))*dz1w(k)
        flv(i,j,k) = 0.5*(Flux(i,j,k+1,1,3) + Flux(i,j,k,1,3))
      end do
    end do
  end do
!
  if (out_diagu) diag(4)%u = diag(4)%u + cdiag*eh
  if (out_diagu) diag(5)%u = diag(5)%u + cdiag*ev
  if (out_fsgs.and.out_u) turbu_diag%ufsgs = -flv
!
  winds%u(it_start:it_end,jt_start:jt_end,1:nz-1) = winds%u(it_start:it_end,jt_start:jt_end,1:nz-1) + &
  		eh(it_start:it_end,jt_start:jt_end,1:nz-1) + ev(it_start:it_end,jt_start:jt_end,1:nz-1)
!
!  V velocity
!
#ifdef MODEL_3D
  eh = 0.0
  ev = 0.0
  flv = 0.0
  do k=1,nz-1
    do j=jt_start,jt_end
      do i=it_start,it_end
	eh(i,j,k) = (Flux(i+1,j,k,2,1) - Flux(i,j,k,2,1))*dx1
!
	eh(i,j,k) = eh(i,j,k) + (Flux(i,j+1,k,2,2) - Flux(i,j,k,2,2))*dy1
!
	ev(i,j,k) = (Flux(i,j,k+1,2,3) - Flux(i,j,k,2,3))*dz1w(k)
        flv(i,j,k) = 0.5*(Flux(i,j,k+1,2,3) + Flux(i,j,k,2,3))
      end do
    end do
  end do
!
#ifdef CHANNEL
  if (jp_start == 1) then
    eh(:,jt_start,:) = 0.
    ev(:,jt_start,:) = 0.
  endif
#endif
!
  if (out_diagv) diag(4)%v = diag(4)%v + cdiag*eh
  if (out_diagv) diag(5)%v = diag(5)%v + cdiag*ev
  if (out_fsgs.and.out_v) turbu_diag%vfsgs = -flv
!
  winds%v(it_start:it_end,jt_start:jt_end,1:nz-1) = winds%v(it_start:it_end,jt_start:jt_end,1:nz-1) + &
  		eh(it_start:it_end,jt_start:jt_end,1:nz-1) + ev(it_start:it_end,jt_start:jt_end,1:nz-1)
#endif
!
!  W velocity
!
  eh = 0.0
  ev = 0.0
  do k=2,nz-1
    do j=jt_start,jt_end
      do i=it_start,it_end
	eh(i,j,k) = (Flux(i+1,j,k,3,1) - Flux(i,j,k,3,1))*dx1
!
#ifdef MODEL_3D
	eh(i,j,k) = eh(i,j,k) + (Flux(i,j+1,k,3,2) - Flux(i,j,k,3,2))*dy1
#endif
!
	ev(i,j,k) = (Flux(i,j,k+1,3,3) - Flux(i,j,k,3,3))*dz1(k)
      end do
    end do
  end do
!
  if (out_diagw) diag(4)%w = diag(4)%w + cdiag*eh
  if (out_diagw) diag(5)%w = diag(5)%w + cdiag*ev
!
  winds%w(it_start:it_end,jt_start:jt_end,2:nz-1) = winds%w(it_start:it_end,jt_start:jt_end,2:nz-1) + &
  		eh(it_start:it_end,jt_start:jt_end,2:nz-1) + ev(it_start:it_end,jt_start:jt_end,2:nz-1)
!
!  Deallocate
!
  call dealloc ( dij )
  call dealloc ( Flux )
  call dealloc ( eh )
  call dealloc ( ev )
  call dealloc ( flv )
  call dealloc ( fkv )
  call dealloc ( fkh )
!
!----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminate efud')
!
return
end 
!
!   ===================================================
    subroutine eff0 ( dens, statel, hydromtrl )
!   ===================================================
!
  type (atm_state), intent(in) :: statel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
  real, dimension(ip_start:,jp_start:,:), intent(in) :: dens
!
  integer :: i, k, h
  real :: func, es0(nz)
  real, dimension(nz) :: zero
!
  real, allocatable, dimension(:,:) :: flux
  real, allocatable, dimension(:,:,:) :: fkv
!
  if (verbose > 1) call write_debug('Enter eff0')
!
!----------------------------------------------------------!
!                 Retrieve scalar diffusivity              !
!----------------------------------------------------------!
!
!  Allocate
!
  call alloc ( fkv )
  call alloc ( flux )
!
  zero = 0.
!
#ifdef ISENTROPIC
  es0 = pt0
#else
  es0 = mse0
#endif
!
!  If Prandtl lower than 0, diffusion is applied in surface layer only
!
  if ( pran < 0.0 ) then
    do k = 1, nz
      func = max(1. - exp(6.*(z0(k) - zdec)/zdec),0.)
      fkv(:,:,k) = dens(:,:,k)*turbu%fkv(:,:,k)*func/abs(pran) 
    enddo
  else
    fkv = dens*turbu%fkv/abs(pran)
  endif
!
!----------------------------------------------------------!
!                       SGS fluxes                         !
!----------------------------------------------------------!
!
!  For Sensible energy
!
  call eff (1, dens, -surf%esflux, es0, states%es, statel%es, fkv)
!
!  For qt
!
  call eff (2, dens, -surf%qvflux, qt0, states%qt, statel%qt, fkv)
!
!  For hydrometeors
!
  flux = 0.
  if (micro_dif) then
    do h = 1, nhydro
!  
!  h Mixing ratio  
!
      call eff (2+h, dens, flux, zero, hydromtrs(h)%q, hydromtrl(h)%q, fkv)
!  
!  h number concentration  
!
#ifdef SEIFERT
      if (moments == 2) call eff (2+nhydro+h, dens, flux, zero, hydromtrs(h)%n, hydromtrl(h)%n, fkv)
#endif
!  
!  h rimed liquid mass  
!
#ifdef SEIFERT
      if ( lmicro > 3 .and. h > ice ) call eff (2+nhydro+h, dens, flux, zero, hydromtrs(h)%w, hydromtrl(h)%w, fkv)
#endif
    enddo 
  endif
!
!  Deallocate
!
  call dealloc ( fkv )
  call dealloc ( flux )
!
!----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminate eff0')
!
return
end
!
! ====================================
  Subroutine effs ( flag, dens, x, es, flux_s )
! ====================================
!
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: dens, es, x
  real, dimension(ip_start:,jp_start:),  intent(IN) :: flux_s
!
  integer :: k, flag
  real :: func
  real, dimension(nz) :: zero
  real, allocatable, dimension(:,:,:) :: fkv, ec
!
!================================================
!
!  Calculation of subgrid scale fluxes for any additional
!  prognostic scalar (chemicals, aerosols, ice nuclei).
!
!================================================
!
!  Allocate
!
  call alloc ( fkv )
  call alloc ( ec )
!
  zero = 0.0
!
!  Get turublent viscosity
!
  if ( pran < 0.0 ) then
    do k = 1, nz
      func = max(1. - exp(6.*(z0(k) - zdec)/zdec),0.)
!
      fkv(:,:,k) = dens(:,:,k)*turbu%fkv(:,:,k)*func/abs(pran) 
    enddo
  else
    fkv = dens*turbu%fkv/abs(pran)
  endif
!
  call eff ( flag, dens, flux_s, zero, ec, x, fkv )
!
  es = es + ec
!
!  Deallocate
!
  call dealloc ( fkv )
  call dealloc ( ec )
!
!----------------------------------------------------------!
!	
return
end
!
!   ================================================
    subroutine eff ( flag, dens, flux_s, x0, ec, x, nut )
!   ================================================

! ----------------------------------------------
! --- Subroutine calculating eddy terms of 
! ---            scalar variables
! ----------------------------------------------
!
  integer :: i, j, k, i1, i2, j1, j2, flag
  real  :: dx1, dy1
  real, dimension(nz) :: f_an
!
  real, dimension(:), intent(in) :: x0
  real, dimension(ip_start:,jp_start:), intent(in) :: flux_s
  real, dimension(ip_start:,jp_start:,:), intent(in) :: dens, x, nut
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: ec
!
  real, dimension(:,:,:), allocatable :: xs, ez, flz, nuh
  real, dimension(:,:,:), allocatable :: gradx, grady, gradz
!
!----------------------------------------------------------!
!                     Initialisations                      !
!----------------------------------------------------------!
!
!  Allocate
!
  call alloc ( ez )
  call alloc ( xs )
  call alloc ( gradx )
  call alloc ( grady )
  call alloc ( gradz )
  call alloc ( flz )
  call alloc ( nuh )
!
  dx1  = 1./dx
#ifdef MODEL_3D
  dy1  = 1./dy
#endif
!
  if (diff /= 0.) then
    do k = 1, nz
      xs(:,:,k) = x(:,:,k) - x0(k)
    enddo
  else
    xs = x
  endif
!
!  Anisotropy factor
!
  f_an = 1.
  if (anis_k) then
    do k = 1, nz
      f_an(k) = ( (dx*dy)**(1./2.) / (dz*fdz0(k)) )**2.
    enddo
  endif
!
  do k = 1, nz
    nuh(:,:,k) = f_an(k)*nut(:,:,k)
  enddo
!
!----------------------------------------------------------!
!          	  Calculate scalar gradient                !
!----------------------------------------------------------!
!
  if (diff_ord == 2) then
    call grad_quick (xs, gradx, grady, gradz)
  else
    call grad_cus (xs, gradx, grady, gradz)
  endif
!
!  Diagnostics
!
  if (out_diagt .and. out_grad .and. flag==1) diag(1)%grad%pt = gradx
#ifdef MODEL_3D
  if (out_diagt .and. out_grad .and. flag==1) diag(2)%grad%pt = grady
#endif
  if (out_diagt .and. out_grad .and. flag==1) diag(3)%grad%pt = gradz
!
  if (out_diagq .and. out_grad .and. flag==2) diag(1)%grad%qt = gradx
#ifdef MODEL_3D
  if (out_diagq .and. out_grad .and. flag==2) diag(2)%grad%qt = grady
#endif
  if (out_diagq .and. out_grad .and. flag==2) diag(3)%grad%qt = gradz
!
!----------------------------------------------------------!
!              Calculate scalar SGS fluxes                 !
!----------------------------------------------------------!
!
  i1=it_start-1
  i2=it_end
  j1=1
  j2=1
#ifdef MODEL_3D
  j1=jt_start-1
  j2=jt_end
#endif
!
!  Fluxes (stored in gradient)
!
  do k = 1, nz-1
    do j = j1, j2
      do i = i1, i2
        gradx(i,j,k) = 0.5*(nuh(i+1,j,k) + nuh(i,j,k)) * gradx(i+1,j,k)
!
#ifdef MODEL_3D
        grady(i,j,k) = 0.5*(nuh(i,j+1,k) + nuh(i,j,k)) * grady(i,j+1,k)
#endif 
!
        gradz(i,j,k) = 0.5*(nut(i,j,k+1) + nut(i,j,k)) * gradz(i,j,k+1)
      enddo
    enddo
  enddo
!
!  In X & Y
!
  do k = 1, nz-1
    do j = jt_start,jt_end
      do i = it_start,it_end
        ec(i,j,k) = (gradx(i,j,k) - gradx(i-1,j,k))*dx1
!
#ifdef MODEL_3D
        ec(i,j,k) = ec(i,j,k) + (grady(i,j,k) - grady(i,j-1,k))*dy1
#endif 
      enddo
    enddo
  enddo
!
!  In Z
!
  do k = 2, nz-1
    do j = jt_start,jt_end
      do i = it_start,it_end
        ez(i,j,k) = (gradz(i,j,k) - gradz(i,j,k-1))*dz1w(k)
	flz(i,j,k) = 0.5*(gradz(i,j,k) + gradz(i,j,k-1))
      enddo
    enddo
  enddo
!
  ez(:,:,1) = (gradz(:,:,1) - flux_s)*dz1w(1)
  flz(:,:,1) = 0.5*(gradz(:,:,1) + flux_s)
!
!----------------------------------------------------------!
!          	         Assemble    		           !
!----------------------------------------------------------!
!
  ec = ec + ez
!
!  Non-conservative case: divide by density afterwards
!
#ifndef CONSERVATIVE
  ec = ec / dens
#endif
!
!  Diagnostics
!
  if (out_diagt .and. flag==1) diag(5)%pt = diag(5)%pt + cdiag*ec 
  if (out_diagq .and. flag==2) diag(5)%qt = diag(5)%qt + cdiag*ec 
  if (out_diagl .and. flag==3) diag(5)%qc = diag(5)%qc + cdiag*ec 
  if (out_diagr .and. flag==4) diag(5)%qr = diag(5)%qr + cdiag*ec 
  if (out_diagi .and. flag==5) diag(5)%qi = diag(5)%qi + cdiag*ec 
!
  if (out_diags .and. flag > 2+2*nhydro) then
    diag(5)%sca(:,:,:,flag-2*(nhydro+1)) = diag(5)%sca(:,:,:,flag-2*(nhydro+1)) + cdiag*ec
  endif
!
  if (out_diaga .and. flag>100 .and. flag<200) diag(4)%aero(flag-100)%n = diag(4)%aero(flag-100)%n + cdiag*ec
  if (out_diaga .and. flag>200) diag(4)%aero(flag-200)%m = diag(4)%aero(flag-200)%m + cdiag*ec
!
  if (out_fsgs) then
    if (out_pt .and. flag==1) turbu_diag%ptfsgs = -flz 
    if (out_qt .and. flag==2) turbu_diag%qtfsgs = -flz 
  endif
!
!  Deallocate
!
  call dealloc ( ez )
  call dealloc ( xs )
  call dealloc ( gradx )
  call dealloc ( grady )
  call dealloc ( gradz )
  call dealloc ( flz )
  call dealloc ( nuh )
!
!----------------------------------------------------------!
!
return                                                                   
end
!
!   ================================================
    subroutine calc_smago ( N2, S2 )
!   ================================================

! ----------------------------------------------
! --- Subroutine calculating eddy viscosity following
!	Smagorinsky-Lilly (including stability correction)
! ----------------------------------------------
!
  integer :: i, j, k
  real  :: zz, Ri, fri
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: N2, S2, delta
!
  character(len=6), parameter :: correction='Lilly' !'Brown'
!
  if (verbose > 1) call write_debug('Enter calc_smago')
!
!----------------------------------------------------------!
!		Characteristic length scale                !
!----------------------------------------------------------!
!
  turbu%fkv  = 0.0
!
  call length_scale ( turbu%ksgs, N2, delta ) 
!
!----------------------------------------------------------!
!          Eddy viscosity using Smagorinsky-Lilly          !
!----------------------------------------------------------!
!
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end
!
!  calculate first buoyancy correction
!	
	Ri = N2(i,j,k) / max(S2(i,j,k), 1.e-12) / abs(pran) 
!
	fri = 1.
	select case (trim(correction))
	case("Lilly") 
	  if (Ri <= Ric) then
	    fri=sqrt(1. - max(Ri,-0.25)/Ric)
	  else
	    fri=0.
	  endif
        case("Brown")
	  if (Ri <= 0.) then
	    fri=sqrt(1. - 4.*max(Ri,-0.25)/Ric)
	  else
	    fri=(1. - min(Ri,Ric)/Ric)**4.
	  endif
	end select
!
!  Eddy viscosity
!
	turbu%fkv(i,j,k) = (C_s*delta(i,j,k))**2. * sqrt( S2(i,j,k) ) * max(fri,0.)
      end do
    end do
  end do
!
!  Boundary conditions
!
  call statebc( turbu%fkv )
!
!  Extra diagnostics
!
  turbu%ksgs = (turbu%fkv / (C_k*delta))**2.
!
!----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminate calc_smago')
!
return                                                                   
end
!
!	==================================================
	subroutine calc_tke ( N2 )
!	==================================================
!
  integer :: i, j, k
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: N2
  real, dimension(:,:,:), allocatable :: delta
!
  if (verbose > 1) call write_debug('Enter calc_tke') 
!
!----------------------------------------------------------!
!                    Set length scale                      !
!----------------------------------------------------------!
!
  call alloc (delta)
!
  call length_scale ( turbu%ksgs, N2, delta ) 
!
!----------------------------------------------------------!
!		  Find eddy viscosities                    !
!----------------------------------------------------------!
!
  turbu%fkv  = 0.0
!   
  turbu%fkv = C_k*delta*sqrt(turbu%ksgs)
!
!  Boundary conditions
!
  call statebc( turbu%fkv )
!
!  Deallocate
!
  call dealloc (delta)
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminate calc_tke')
!
return
end
!
!   ================================================
    subroutine stability ( pressurel, statel, hydromtrl, dij, N2, S2 )
!   ================================================

! ----------------------------------------------
! --- Subroutine calculating stability parameters used
!	in SGS turbulence closure
! ----------------------------------------------
!
  type (atm_state), intent(in) :: statel
  type (atm_pressure), intent(in) :: pressurel
  type (hydrometeor), dimension(1:nhydro), intent(in) :: hydromtrl
  real, dimension(ip_start:,jp_start:,:,:,:), intent(in) :: dij
!
  real, dimension(ip_start:,jp_start:,:), optional, intent(inout) :: N2, S2
!
  integer :: i, j, k, h
  real :: A, qv, tem, sumh, n2l
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: ptv
!
!----------------------------------------------------------!
!               Calculate parameters for K                 !
!----------------------------------------------------------!
!
!  Calculate N2
!
  if (present(N2)) then
!
  call get_ptv ( pressurel, statel, hydromtrl, ptv )
!
  do k=2,nz-1
    do j=jt_start,jt_end
      do i=it_start,it_end
        sumh = 0.
        do h = 1, nhydro
	  sumh = sumh + hydromtrl(h)%q(i,j,k)
	enddo
        qv = statel%qt(i,j,k) - sumh
	tem = thermo%T(i,j,k)
!
!  Need to distinguish between saturated and unsaturated layers (Klemp and Wilhelmson, 1978)
!
	if (lmicro > 0 .and. with_mic) then
	  if ( hydromtr2(drop)%q(i,j,k)+hydromtr2(rain)%q(i,j,k) < 1.e-6 ) then
	    n2l = g / (Pt0(k)*(1. + epsm*qv0(k))) * (ptv(i,j,k+1) - ptv(i,j,k-1)) / (dz*(0.5*fdz0(k+1)+fdz0(k)+0.5*fdz0(k-1)))
	  else
	    A  = (1. + 1.61*epsm*cal_flv(tem)*qv/((cp_a-cv_a)*tem)) / (1. + epsm*cal_flv(tem)**2.*qv/(cp_a*(cp_a-cv_a)*tem**2.))
!
	    n2l = g * (A/(Pt0(k)*(1.+epsm*qv0(k)))*(ptv(i,j,k+1) - ptv(i,j,k-1)) - 	&
	         (hydromtr2(drop)%q(i,j,k+1) - hydromtr2(drop)%q(i,j,k-1) + hydromtr2(rain)%q(i,j,k+1) -	&
	          hydromtr2(rain)%q(i,j,k-1))) / (dz*(0.5*fdz0(k+1)+fdz0(k)+0.5*fdz0(k-1)))
	  endif
	else
	  n2l = g / (Pt0(k)*(1. + epsm*qv0(k))) * (ptv(i,j,k+1) - ptv(i,j,k-1)) / (dz*(0.5*fdz0(k+1)+fdz0(k)+0.5*fdz0(k-1)))	
	endif
!
        N2(i,j,k) = n2l
      end do
    end do
  end do
!
!  Surface and top values
!
  do j=jt_start,jt_end
    do i=it_start,it_end
      N2(i,j,nz) = 0.5*g/(Pt0(nz)*(1. + epsm*qv0(nz))) * (ptv(i,j,nz) - ptv(i,j,nz-1)) / (0.5*dz*(fdz0(nz)+fdz0(nz-1)))
      N2(i,j,1) = 0.5*g/(Pt0(1)*(1. + epsm*qv0(1))) * (ptv(i,j,2) - ptv(i,j,1)) / (0.5*dz*(fdz0(2)+fdz0(1)))
    enddo
  enddo
!
  call statebc ( N2 )
!
  endif
!
!  Calculate S2
!
  if (present(S2)) then
!
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end
	S2(i,j,k) = max(2.*sum(dij(i,j,k,:,:)*dij(i,j,k,:,:)),1.e-12)
      end do
    end do
  end do
!
  call statebc ( S2 )
!
  endif
!
!----------------------------------------------------------!
!
return                                                                   
end
!
! ==================================================
  subroutine tkecomp ( pressurel, windl, statel, hydromtrl, N2, S2 ) 
! ==================================================
!
  type (atm_pressure), intent(in) :: pressurel
  type (atm_winds), intent(in) :: windl
  type (atm_state), intent(in) :: statel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
  real, dimension(ip_start:,jp_start:,:), intent(in) :: N2, S2
!
  integer :: i, j, k
  real :: dt1, l_d, fk, Ri, fri
!
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: zero2d
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: tke, tkes, delta, fkv
!
  if (verbose > 1) call write_debug('Starting tkecomp')
!
!----------------------------------------------------------!
!                    Set parameters                        !
!----------------------------------------------------------!
!
  dt1 = 1./dt0
!
  zero2d = 0.0
  delta = 0.0
  tkes = 0.0
  fkv = 0.0
  tke = turbu%ksgs
!
!  Turbulent scale
!
  call length_scale ( tke, N2, delta ) 
!
!  Turbulent viscosity
!
  fkv = C_k*delta*sqrt(tke) 
!
!----------------------------------------------------------!
!                         SGS Fluxes                       !
!----------------------------------------------------------!
!
  if (with_dif) call eff ( 0, pressurel%dens, zero2d, k0, tkes, tke, 2.*fkv )
!
!----------------------------------------------------------!
!                       Advection                          !
!----------------------------------------------------------!
!
  if (with_adv) then
    call advection_x ( 'k ', dt0, tke, turbu%ksgs, tkes, wind_adv%u, .true., adv=.true. )
!
#ifdef MODEL_3D
    call advection_y ( 'k ', dt0, tke, turbu%ksgs, tkes, wind_adv%v, .true., adv=.true. )
#endif
!
    call advection_z ( 'k ', dt0, tke, turbu%ksgs, tkes, wind_adv%w, .true., adv=.true. )
!
!  Assemble advection tendency
!
#ifdef ADV_SPLIT
    tkes = tkes + (tke - turbu%ksgs)*dt1
#endif
  endif
!
!----------------------------------------------------------!
!                    Sources and sinks                     !
!----------------------------------------------------------!
!
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end
	Ri = N2(i,j,k)/max(S2(i,j,k),1.e-12)
	if (Ri <= Ric) then
	  fri = 1. - max(Ri/abs(pran),-1.)
	else
	  fri = 0.
	endif
!
        tkes(i,j,k) = tkes(i,j,k) + fri*fkv(i,j,k)*S2(i,j,k) - C_e*turbu%ksgs(i,j,k)**(1.5)/delta(i,j,k)
      end do
    end do
  end do
!
!----------------------------------------------------------!
!                Integrate and terminate                   !
!----------------------------------------------------------!
!
!  Integrate
!
  turbu%ksgs = max(turbu%ksgs + dt0*tkes,0.0)
!
!  Boundary conditions      
!
  call statebc( turbu%ksgs )
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminating tkecomp')
!
return
end
!
! ==================================================
  subroutine non_linear ( pressurel, windl, statel, hydromtrl, dij, fkv ) 
! ==================================================
!
  type (atm_pressure), intent(in) :: pressurel
  type (atm_winds), intent(in) :: windl
  type (atm_state), intent(in) :: statel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
!
  real, dimension(ip_start:,jp_start:,:,:,:), intent(inout) :: dij
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: fkv
!
  integer :: i, j, k, l, m
  real :: NL
  real, dimension(1:3,1:3) :: del
  real, dimension(:,:,:), allocatable :: S2, N2, delta
  real, dimension(:,:,:,:,:), allocatable :: rij
!
!----------------------------------------------------------!
!                    Set parameters                        !
!----------------------------------------------------------!
!
!  Allocate and initialise
!
  call alloc ( rij, 3, 3 )
  call alloc ( S2 )
  call alloc ( N2 )
  call alloc ( delta )
!
  del = 0.
  do l = 1, 3
    del(l,l) = 1.
  enddo
!
  call stability ( pressurel, statel, hydromtrl, dij, N2=N2, S2=S2 )
!
  call length_scale ( turbu%ksgs, N2, delta )
!
!----------------------------------------------------------!
!                  Non-linear contribution                 !
!----------------------------------------------------------!
!
  call cal_dij ( windl, rij=rij ) 
!
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
!
        do l = 1, 3
          do m = 1, 3
            NL = (sum(dij(i,j,k,l,:)*dij(i,j,k,:,m)) - 1./6.*S2(i,j,k)*del(l,m)) + (sum(dij(i,j,k,l,:)*rij(i,j,k,:,m)) - sum(rij(i,j,k,l,:)*dij(i,j,k,:,m)))
            dij(i,j,k,l,m) = 2.*dij(i,j,k,l,m) + C_1*NL/sqrt(max(S2(i,j,k),1.e-12))
          enddo
        enddo
!
      enddo
    enddo
  enddo
!
#ifdef TKE
  fkv = C_k*pressurel%dens*delta*sqrt(turbu%ksgs) 
#else
  fkv = (C_s*delta)**2.*pressurel%dens*sqrt(S2)
#endif
!
!  Deallocate
!
  call dealloc ( rij )
  call dealloc ( S2 )
  call dealloc ( N2 )
  call dealloc ( delta )
!
!-----------------------------------------------------------!
!
return
end
!
! ==================================================
  subroutine length_scale ( ksgs, N2, delta ) 
! ==================================================
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: N2, ksgs 
  real, dimension(ip_start:,jp_start:,:),  intent(inout) :: delta
!
  integer :: i, j, k
  real :: l_d, l_z
  real, dimension(1:nz) :: deltas
!
!----------------------------------------------------------!
!                    Set parameters                        !
!----------------------------------------------------------!
!
  deltas = 0.
  delta = 0.
!
  if (.not.anis_k) then
    do k = 1, nz 
#ifdef MODEL_3D
      deltas(k) = (dx*dy*dz*fdz0(k))**(1./3.) 
#else
      deltas(k) = (dx*dz*fdz0(k))**(1./2.)
#endif
    enddo
  else
    deltas = dz*fdz0
  endif
!
!  Turbulent scale limited by Deardorff scale and surface
!
  do k = 1, nz
    do j = jp_start, jp_end
      do i = ip_start, ip_end
        l_z = vk*(z0(k) + zrough) 
        l_d = 1e8
#ifdef TKE
        if (N2(i,j,k) > 1.e-12) l_d = C_n*sqrt(max(ksgs(i,j,k),1.e-12) / N2(i,j,k))
#endif
!
        delta(i,j,k) = 1. / sqrt(1./deltas(k)**2. + 1./l_d**2. + 1./l_z**2.)
      enddo
    enddo
  enddo
!
!-----------------------------------------------------------!
!
return
end
!
end module subgrid

