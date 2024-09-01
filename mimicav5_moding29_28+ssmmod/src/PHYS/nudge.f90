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
!  nudge.f90                 
!
!  Purpose:
!	Nudging of main transported quantities	
!
!  Author
!	Julien Savre
!	MISU
!
! ================================================================

!	==================================================
	subroutine nudge_uv
!	==================================================
!
! ---------------------------------
! --- Calculate nudging for horizontal wind components
! ---------------------------------
!
USE gridno
USE shared_data	
USE shared_wind
USE shared_pressure
USE shared_diag
USE averages
!
IMPLICIT NONE
!
  integer :: i, j, k
  real  :: zz
  real, dimension(nz)  :: cuv, uav, vav,un,vn
!
! ---------------------------------------------------------!
!
!  Select case for the expression of nudging coefficients
!
  cuv = 0.
  un = u00
  vn = v00
!
  select case (casename)
!
    case ('ISDAC')
      do k = 1, nz
        if (z0(k) <= 825.) then
          cuv(k) = 1./tnudg * 0.5*(1. - cos(pi*z0(k)/825.))
        else
          cuv(k) = 1./tnudg
        endif
      enddo    
!
    case ('MOCCHA')
      do k = 1, nz
        if (z0(k) <= 1600.) then
          cuv(k) = 1./tnudg * 0.5*(1. - cos(pi*z0(k)/1600.))
        else
          cuv(k) = 1./tnudg
        endif
      enddo    
!
    case('MIDLAT')
      do k = 1, nz
        cuv(k) = 0.5*(1. + tanh((50000. - p0(k)) / 22000.)) / tnudg
      enddo
!
    case ('RCE2','RCE4','RCE8','RCE12','RCEMIP','TRANSIT','RCERH')
      cuv = 1./tnudg
!      
    case ('SE_ATLANTIC')
      un = u_force     
      vn = v_force  
      cuv = coef_nudge2 !inverse of the nudging time scale
  end select
!
!  Calculate nudging terms
!
  call horav ( wind%u, uav )
#ifdef MODEL_3D
  call horav ( wind%v, vav )
#endif
!
  do k = 1, nz
    winds%u(:,:,k) = winds%u(:,:,k) - cuv(k)*(uav(k) - den0(k)*un(k))
#ifdef MODEL_3D
    winds%v(:,:,k) = winds%v(:,:,k) - cuv(k)*(vav(k) - den0(k)*vn(k))
#endif
  enddo        
!
  if (out_diagu) then
    do k = 1, nz
      diag(8)%u(:,:,k) = diag(8)%u(:,:,k) - cdiag*cuv(k)*(uav(k) - den0(k)*un(k))
    enddo
  endif
!
#ifdef MODEL_3D
  if (out_diagv) then
    do k = 1, nz
      diag(8)%v(:,:,k) = diag(8)%v(:,:,k) - cdiag*cuv(k)*(vav(k) - den0(k)*vn(k))
    enddo
  endif
#endif
!
! ---------------------------------------------------------!
!
return
end  
!
!	==================================================
	subroutine coriolis
!	==================================================
!
! ---------------------------------
! --- Calculate nudging for horizontal wind components
! ---------------------------------
!
USE gridno
USE averages
USE shared_data	
USE shared_wind
USE shared_pressure
USE shared_diag
!
IMPLICIT NONE
!
  integer :: i, j, k
  real, dimension(1:nz) :: ugg, vgg, uav, vav
  real  :: l, windm
  real, dimension(jp_start:jp_end) :: f
!
! ---------------------------------------------------------!
!
#if (defined MODEL_3D)
!
#ifdef CHANNEL
  do j = jt_start, jt_end
    l = abs(real(j-3) - 0.5*real(ny-5))*dx/111100.
    f(j) = f2omega*sin(pi*(ctr_lat+l)/180.)
  enddo
#else
  f = f2omega*sin(pi*ctr_lat/180.)
#endif
!
  select case (trim(casename))
!  
  case ('BOMEX')
!   
    vgg = 0. 
    do k = 1, nz
      ugg(k) = -10. + 1.8e-3 * z0(k)
    enddo
!
    call horav(wind%u, uav)
    call horav(wind%v, vav)
!
    do k=1,nz
      winds%u(:,:,k) = winds%u(:,:,k) + 0.376e-4*(vav(k) - den0(k)*vgg(k))
      winds%v(:,:,k) = winds%v(:,:,k) - 0.376e-4*(uav(k) - den0(k)*ugg(k))
!
      if (out_diagu) diag(6)%u(:,:,k) = diag(6)%u(:,:,k) + cdiag*0.376e-4*(vav(k) - den0(k)*vgg(k))
      if (out_diagv) diag(6)%v(:,:,k) = diag(6)%v(:,:,k) - cdiag*0.376e-4*(uav(k) - den0(k)*ugg(k))
    enddo
! 
  case default
!    
    ugg = u00
    vgg = v00
!
    do k=1,nz
      do j=jt_start,jt_end
        do i=it_start,it_end
          windm = 0.25*(wind%v(i,j,k) + wind%v(i,j+1,k) + wind%v(i-1,j,k) + wind%v(i-1,j+1,k))
	  winds%u(i,j,k) = winds%u(i,j,k) + f(j)*(windm - den0(k)*vgg(k))
          if (out_diagu) diag(6)%u(i,j,k) = diag(6)%u(i,j,k) + cdiag*f(j)*(windm - den0(k)*vgg(k))
!
          windm = 0.25*(wind%u(i,j,k) + wind%u(i+1,j,k) + wind%u(i,j-1,k) + wind%u(i-1,j+1,k))
	  winds%v(i,j,k) = winds%v(i,j,k) - f(j)*(windm - den0(k)*ugg(k))
          if (out_diagv) diag(6)%v(i,j,k) = diag(6)%v(i,j,k) - cdiag*f(j)*(windm - den0(k)*ugg(k))
        enddo
      enddo
    enddo
!
  end select
!
#endif
!
! ---------------------------------------------------------!
!
return
end  
!
!	==================================================
	subroutine nudge_x ( c3, x0, x, xs )
!	==================================================
!
! ---------------------------------
! --- Calculate nudging for horizontal wind components
! ---------------------------------
!
USE gridno
USE shared_pressure
USE shared_data	
USE shared_diag
USE averages
!
IMPLICIT NONE
!
  integer :: k
  character(len=3) :: c3
  real, dimension (ip_start:ip_end,jp_start:jp_end,1:nz) :: xs, x, nud
!
  real, dimension(nz)  :: x0, xn, cxx, xav
  real  :: tnudge, zz, ddz=200.
!
! ---------------------------------------------------------!
!
!  Select case for the expression of nudging coefficients
!
  cxx = 0.
  xn = x0  !nudging profile (from the forcing file)
  nud(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.
!
  select case (trim(casename))
!
    case ('ISDAC')
      do k = 1, nz
        if ( z0(k) < 1200. ) then
          cxx(k) = 0.
	else if ( z0(k) >= 1200. .and. z0(k) <= 1500. ) then
	  cxx(k) = 2./tnudg * 0.5*( 1. - cos(pi*(z0(k)-1200.)/300.) )
	else
	  cxx(k) = 2./tnudg
	endif
      enddo    
!
    case('MIDLAT')
      do k = 1, nz
        cxx(k) = 0.5*(1. + tanh((50000. - p0(k)) / 22000.)) / tnudg
      enddo
!     
!    case ('RCEMIP','COLUMN')   
!      if (c3=='sca' .and. maxval(s_scal0) == 0.) then
!        do k = 1, nz
!          if ( z0(k) >= 0. .and. z0(k) < 15000. ) then
!            tnudge = 86400.*(5. + 95.*z0(k)/15000.)
!            cxx(k) = 1./tnudge
!          else if ( z0(k) >= 15000. ) then
!            tnudge = 86400.*max(100. - 99.9*(z0(k) - 15000.)/30000., 0.1)
!            cxx(k) = 1./tnudge
!          endif
!        enddo 
!      endif
!
    case ('RCERH')
      if (c3 == 'qt') then
        cxx = 3./tnudg
      endif   
! 
    case ('SE_ATLANTIC')
      cxx = coef_nudge1 !nudging time scale
      if (c3=='pt ') xn = ta_force
      if (c3=='qt ') xn = qt_force  
      FORALL(k=1:nz) xn(k) = xn(k) * den0(k) !multiplying nudging profile (qt or PT) by the density

  end select
!
!  Calculate nudging terms
!
  call horav ( x, xav )
!
  do k = 1, nz
    nud(:,:,k) = cxx(k)*(xn(k) - xav(k))
  enddo 
!
!  Terminate
!
  xs(ip_start:ip_end,jp_start:jp_end,1:nz) = xs(ip_start:ip_end,jp_start:jp_end,1:nz) +   &
  	nud(ip_start:ip_end,jp_start:jp_end,1:nz)
!
  if (out_diagt .and. c3=='pt ') diag(8)%pt = diag(8)%pt + cdiag*nud
  if (out_diagq .and. c3=='qt ') diag(8)%qt = diag(8)%qt + cdiag*nud
  if (out_diags .and. c3=='sca') diag(8)%sca(:,:,:,1) = diag(8)%sca(:,:,:,1) + cdiag*nud
!
! ---------------------------------------------------------!
!
return
end
!
! ==================================================
  subroutine ls_advec ( flag, x, s )
! ==================================================
!
! ---------------------------------
! --- Calculate large scale scalar advection
! ---------------------------------
!
USE gridno
USE shared_data	
USE shared_pressure
USE shared_thermo
USE shared_diag
USE averages
!
IMPLICIT NONE
!
  integer :: i1, i2, j1, j2, k1, k2
  integer :: i, j, k, flag, it
!
  real, dimension (ip_start:ip_end,jp_start:jp_end,1:nz) :: s, x
  real, dimension(nz)  :: lsadv, pav
  real  :: pp, lsa1, lsa2
!
! ---------------------------------------------------------!
!
  it = time/10800.
!
!  Select case: First RICO
!
  if ( casename(1:4) == 'RICO' ) then
!
    select case (flag)
!
!  Large scale advection of temperature
!
      case (1)
        lsadv = -2.5 / 86400.
!     
!  Large scale advection of moisture
! 
      case (2)
        do k = 1, nz
          if ( z0(k) < 2980. ) then
	    lsadv(k) = -(1.e-3 - 1.3456e-3 * z0(k)/2980.) / 86400.
	  else
  	    lsadv(k) = 0.3456e-3 / 86400.  
	  endif
	enddo
!      
    end select
!
!  Select case: Second BOMEX
!
  else if ( trim(casename) == 'BOMEX' ) then
!
    select case (flag)
!
!  Large scale advection of temperature
!
      case (1)
        do k = 1, nz
          if ( z0(k) < 1500. ) then
	    lsadv(k) = -2. / 86400.
	  else
  	    lsadv(k) = -2./86400. + 2./86400. * (z0(k) - 1500.) / 1500.
	  endif
	enddo
!     
!  Large scale advection of moisture
! 
      case (2)
        do k = 1, nz
          if ( z0(k) < 300. ) then
	    lsadv(k) = -1.2e-8
	  else if ( z0(k) <= 500. ) then
  	    lsadv(k) = -1.2e-8 + 1.2e-8 * (z0(k) - 300.) / 200. 
	  else
	    lsadv(k) = 0.
	  endif
	enddo
!      
    end select
!
!  Default case
!
  else
!
    select case (flag)
!
!  Large scale advection of temperature
!
      case (1)
        call horav (pressure%p, pav)
        do k = 1, nz
          pp = 1.e+5*pav(k)**(cp/Ra)
	  do i = 1, 7
  	    if (0.01*pp <= lsp(i) .and. 0.01*pp > lsp(i+1)) then
              lsa1 = lstem0(i) + (time - real(it)*10800.) / 10800. * (lstem(i) - lstem0(i))	
              lsa2 = lstem0(i+1) + (time - real(it)*10800.) / 10800. * (lstem(i+1) - lstem0(i+1))
              lsadv(k) = lsa1 + (0.01*pp - lsp(i)) / (lsp(i+1) - lsp(i)) * (lsa2 - lsa1)
	    endif
	  enddo 
        enddo    
        lsadv = lsadv * pav
!     
!  Large scale advection of moisture
! 
      case (2)
        call horav (pressure%p, pav)
        do k = 1, nz
          pp = 1.e+5*pav(k)**(cp/Ra)
	  do i = 1, 7
  	    if (0.01*pp <= lsp(i) .and. 0.01*pp > lsp(i+1)) then
              lsa1 = lsqt0(i) + (time - real(it)*10800.) / 10800. * (lsqt(i) - lsqt0(i))
	      lsa2 = lsqt0(i+1) + (time - real(it)*10800.) / 10800. * (lsqt(i+1) - lsqt0(i+1))
              lsadv(k) = lsa1 + (0.01*pp - lsp(i)) / (lsp(i+1) - lsp(i)) * (lsa2 - lsa1)
	    endif
	  enddo 
        enddo          
!      
    end select
!
  endif
!
!  Calculate lsadv
!
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        s(i,j,k) = s(i,j,k) + lsadv(k)
!
        select case (flag)
          case (1)
            if (out_diagt) diag(8)%pt(i,j,k) = diag(8)%pt(i,j,k) + cdiag*lsadv(k)
          case (2)
 	    if (out_diagq) diag(8)%qt(i,j,k) = diag(8)%qt(i,j,k) + cdiag*lsadv(k)
        end select
      enddo
    enddo
  enddo
!
! ---------------------------------------------------------!
!
return
end 

