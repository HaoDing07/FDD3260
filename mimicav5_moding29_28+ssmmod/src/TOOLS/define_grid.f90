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
!  INITIAL.F                   
!
!  Purpose:
!	Calculation of variable cell sizes in vertical direction			  
!
!  Author
!	Julien Savre, MISU
!
! ================================================================

	subroutine define_grid (nnz, dzz, fdz, zz)

!    =================================================

USE gridno
USE shared_data
!
IMPLICIT NONE
!
  integer  :: i, k, io, nnz
  real :: dzz
  real, dimension(1:nnz) :: fdz, zz
  character(len=15) :: refine_name
  logical :: ex=.false., op
!                                                          !
!----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Starting define_grid')
!
!  Grid file present?
!
  inquire ( FILE=trim(gridfile), EXIST=ex )
  if (ex) goto 101
!
!  Else
!
#ifndef FINE
  fdz(:) = 1.
  goto 100
#endif
!
  select case (trim(casename))
    case ('DYCOMS')
      refine_name = 'DYCOMS'
    case ('SE_ATLANTIC')
      refine_name = 'SE_ATLANTIC'
    case ('RCE2','RCE4','RCE8','RCE12','TRANSIT','MIDLAT','CP_IDEAL','ISDAC','MOCCHA')
      refine_name = 'BL'
    case ('ASCOS')
      refine_name = 'sin2'
    case ('GERMANY')
      refine_name = 'geo'
    case('SPIRAL','CONE')
      refine_name = 'cos'
    case default
      refine_name = 'default'
    end select
!
!----------------------------------------------------------!
!                Calculate cell size ratio                 ! 
!----------------------------------------------------------!
!
  select case (trim(refine_name))
!
    case ('constant')
      call cst_refine (nnz, dzz, fdz)
!
    case ('sin2')
      call sin2_refine (nnz, dzz, fdz)
!
    case ('cos')
      call cos_refine(nnz, dzz, fdz)
!
    case ('BL')
      call bl_refine (nnz, dzz, fdz)
!
    case ('DYCOMS')
      !call DYCOMS_refine (nnz, dzz, fdz)
      call DYCOMS_refine_b (nnz, dzz, fdz)  

    case ('SE_ATLANTIC') 
      call SE_Atlantic_refine(nnz,dzz,fdz)  
!
    case ('sin')
      call sin_refine (nnz, dzz, fdz)
!
    case ('geo')
      call geo_refine (nnz, dzz, fdz)
!
    case ('zones')
      call zones_refine (nnz, dzz, fdz)
!
    case default
      write(7,1304) 
1304  format(1x/6x,'Unrecognized refinement method, set constant cell size',/) 
      fdz(:) = 1.
!
  end select
!
100 continue
!
!----------------------------------------------------------!
!           Write vertical dimension output file           ! 
!----------------------------------------------------------!
!
  kbl = 0
  zz(1) = 0.5*dzz*fdz(1)
  do k = 2, nnz
    zz(k) = zz(k-1) + 0.5*dzz*(fdz(k-1)+fdz(k))
    if (kbl == 0 .and. zbl > 0. .and. zz(k) >= zbl) kbl = k
  enddo
  zmax = zz(nnz)
!
!----------------------------------------------------------!
!                Vertical grid from file 	           ! 
!----------------------------------------------------------!
!
101 continue
  if (ex) then 
!
!  Open and read
!
    do i = 10, 100
      inquire ( unit=i, opened=op )
      if (.not.op) exit
    enddo
    
    open ( i, FILE=trim(gridfile), FORM='formatted', STATUS='OLD' )

    kbl = 0  
    
    do k = 1, nnz
      read (i,*,iostat=io) zz(k)
      if (kbl == 0 .and. zbl > 0. .and. zz(k) >= zbl) kbl = k
      if (io /= 0) write(7,1305)
    enddo
    close ( i )
    zmax = zz(nnz)
!
!  Get grid spacing
!
    fdz(1) = 2.*zz(1)/dzz
    do k = 2, nnz
      fdz(k) = 2.*(zz(k) - zz(k-1))/dzz - fdz(k-1)
    enddo
!
  endif
1305  format(1x/6x,'Input grid file not large enough. Complete it or decrease nz.',/) 
!
!----------------------------------------------------------!
!           Calculate the inverse fractions of dz          ! 
!----------------------------------------------------------!
!
  do k = 1, nnz
    dzw(k) = dzz*fdz(k)
    dz1w(k) = 1. / (dzz * fdz(k))
  enddo
!
  do k = 1, nnz-1
    dz1(k) = 1. / (dzz * 0.5*(fdz(k) + fdz(k+1)))
  enddo
  dz1(nnz) = dz1(nnz-1)
!
!  Tropopause level
!
  ktrop = k+1
!                                                          !
!----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminating define_grid')
!
  return
end
!
!    =================================================
	subroutine cst_refine (nnz, dzz, fdz)
!    =================================================
!
USE gridno
USE shared_data
!
IMPLICIT NONE
!
  integer  :: k, kk, nnz
  real, dimension(1:nnz) :: fdz
  real     :: dzz, z, a=1.25
!                                                          !
!----------------------------------------------------------!
!
!  Write informations in output file
!
	write(7,1104) z1, z2
1104	format(8x,'Refined zone boundaries = ',2i5/)
!
	write(7,1204) sratio
1204	format(8x,'Size ratio of the small cells = ',i5/)
!
!----------------------------------------------------------!
!                Calculate cell size ratio                 ! 
!----------------------------------------------------------!
!
!  First zone: from bottom up to z1, dz(k) = dz0
! 
  z = 0.0
  k = 1
  do while (z <= z1)
    fdz(k) = 1.
    z = z + dzz*fdz(k)
    k = k+1
  enddo
!
!  Second zone: Between z1 and z2, dz(k) = dz0*sratio
! 
  do while (z <= z2)
    fdz(k) = sratio
    z = z + dzz*fdz(k)
    k = k+1
  enddo
!
!  Third zone: From z2 up to the top, dz(k) = dz0
! 
  do kk = k, nnz
    fdz(kk) = fdz(kk-1)*a
    z = z + dzz*fdz(kk)
  enddo
!                                                          !
!----------------------------------------------------------!
!
  return
end
!
!    =================================================
	subroutine zones_refine (nnz, dzz, fdz)                                        
!    =================================================
!
USE gridno
USE shared_data
!
IMPLICIT NONE
!
  integer  :: k, nnz
  real, dimension(1:nnz) :: fdz
  real     :: z, dzz
!                                                          !
!----------------------------------------------------------!
!
!  Write informations in output file
!
	write(7,1104) z1, z2
1104	format(8x,'Refined zone boundaries = ',2f7.1/)
!
	write(7,1204) sratio
1204	format(8x,'Size ratio of the small cells = ',f7.2/)
!
!----------------------------------------------------------!
!                Calculate cell size ratio                 ! 
!----------------------------------------------------------!
!
!  First zone: from bottom up to z1, dz(k) = dz0
! 
  z = 0.0
  k = 1
  do while (z <= z1)
    fdz(k) = sratio
    z = z + dzz*fdz(k)
    k = k+1
  enddo
!
!  Second zone: Between z1 and z2, dz(k) = dz0*sratio
! 
  do while (k <= nnz)
    fdz(k) = 1.
    z = z + dzz*fdz(k)
    k = k+1
  enddo
!                                                          !
!----------------------------------------------------------!
!
  return
end
!
!    =================================================
	subroutine bl_refine (nnz, dzz, fdz)                                        
!    =================================================
!
USE gridno
USE shared_data
!
IMPLICIT NONE
!
  integer  :: k, l, m, nnz
  real     :: z, dzm, dzz, S, q, dq
  real, dimension(1:nnz) :: fdz
!                                                          !
!----------------------------------------------------------!
!
!  Write informations in output file
!
	write(7,1104) z1, ztop
1104	format(8x,'Refined boundary layer = ',2f7.1/)
!
	write(7,1204) sratio
1204	format(8x,'Size ratio of the small cells = ',f7.2/)
!
!----------------------------------------------------------!
!                Calculate cell size ratio                 ! 
!----------------------------------------------------------!
!
!  First zone: from bottom up to z1, dz(k) = dz0
! 
  z = 0.0
  k = 1
  do while (z < z1)
    fdz(k) = sratio
    z = z + dzz*fdz(k)
    k = k+1
  enddo
!
!  Second zone: Between z1 and z2, dz(k) = dz0*sratio
!
  S = ztop - z
  m = nnz - k 
  q = 1.2
  dq = 1.
  do while (abs(dq) > 1.e-3)
    dq = (1. - q)*(1. - q**(m+1) - (1. - q)*S/dzz/sratio) / (1 - (m+1)*q**m + m*q**(m+1))
    q = q - dq
  enddo
!
  dzm = sratio*dzz
  do while (k <= nnz)
    fdz(k) = dzm/dzz
    z = z + dzz*fdz(k)
    dzm = q*dzm
    k = k+1
  enddo
!                                                          !
!----------------------------------------------------------!
!
  return
end
!
!    =================================================
	subroutine sin2_refine (nnz, dzz, fdz)
!    =================================================
!
USE gridno
USE shared_data
!
IMPLICIT NONE
!
  integer  :: k, i, nzt, nnz
  real, dimension(1:nz) :: fdz
  real     :: dzz, dzm, a, ierr, q, q0, zzt, z
!                                                          !
!----------------------------------------------------------!
!
!  Write informations in output file
!
	write(7,1104) z1, z2
1104	format(8x,'Refined zone boundaries = ',2f8.3/)
!
	write(7,1204) sratio
1204	format(8x,'Size ratio of the small cells = ',1f8.3/)
!
!----------------------------------------------------------!
!                Calculate cell size ratio                 ! 
!----------------------------------------------------------!
!
  dzm = dzz * sratio
!
!  First zone: from bottom up to nzf1 (height in m), sin2 function
! 
  z = 0.0
  k = 1
  do while (z <= z1)
    fdz(k) = dzm + (dzz - dzm) * sin(pi/z1*z) * sin(pi/z1*z)
    z = z + fdz(k)
    k = k+1
  enddo
!
!  Second zone: Between z1 and z2, dz(k) = cte
! 
  do while (z <= z2)
    fdz(k) = dzm
    z = z + fdz(k)
    k = k+1
  enddo
!
!  Third zone: From z2 up to the top, stretching geometric
! 
  if (k > nnz-8) then
    print*, 'ERROR: the total number of points in the vertical, nz, is not sufficient'
    return
  else
    nzt = nnz - k
    zzt = ztop - (z - dzm)
    a   = dzm / zzt
!
!  Newton to find the geometric factor
!
    q0 = 3.
    ierr = 1.
    do while (ierr > 1.e-4)
      q = q0 - (a*q0**(nzt+1) - q0 + (1.-a)) / (a*(nzt+1)*q0**nzt - 1.)
      ierr = abs(q-q0) 
      q0 = q
    enddo
  endif
!
  do i = k, nnz
    fdz(i) = fdz(i-1) * q
  enddo
  fdz = fdz/dzz
!
  return
  end
!
!    =================================================
	subroutine cos_refine (nnz, dzz, fdz)
!    =================================================
!
USE gridno
USE shared_data
!
IMPLICIT NONE
!
  integer  :: k, nnz
  real     :: dzz
  real, dimension(1:nnz) :: fdz
!                                                          !
!----------------------------------------------------------!
!                Calculate cell size ratio                 ! 
!----------------------------------------------------------!
!
  do k = 1, nnz
    fdz(k) = 0.5*(2. + cos(2.*real(k)/real(nnz)*pi)) 
  enddo
!
  return
  end
!
!    =================================================
	subroutine geo_refine (nnz, dzz, fdz)
!    =================================================
!
USE gridno
USE shared_data
!
IMPLICIT NONE
!
  integer  :: k, nnz
  real, dimension(1:nnz) :: fdz
  real     :: dzz, z, zi
!                                                          !
!----------------------------------------------------------!
!                Calculate cell size ratio                 ! 
!----------------------------------------------------------!
!
  zi = 0.0
  k = 1
  do while (k <= nnz)
    z = dzz*real(k)**sratio
    fdz(k) = z - zi
    zi = z
    k = k+1
  enddo
  fdz = fdz/dz
!
  return
  end
!
!    =================================================
	subroutine sin_refine (nnz, dzz, fdz)
!    =================================================
!
USE gridno
USE shared_data
!
IMPLICIT NONE
!
  integer  :: k, nnz
  real, dimension(1:nnz) :: fdz
  real     :: dzz, b, z
!                                                          !
!----------------------------------------------------------!
!
!  Write informations in output file
!
	write(7,1104) z1, z2
1104	format(8x,'Refined zone boundaries = ',2f8.3/)
!
	write(7,1204) sratio
1204	format(8x,'Size ratio of the small cells = ',1f8.3/)
!
!----------------------------------------------------------!
!                Calculate cell size ratio                 ! 
!----------------------------------------------------------!
!
!  First zone: from bottom up to nzf1 (height in m), sin2 function
! 
  b = 2.5
!
  z = 0.0
  k = 1
  do while (k <= nnz)
    fdz(k) = b*dzz*(1. - 0.7817530784*sin(pi*z/ztop))
    z = z + fdz(k)
    k = k+1
  enddo
  fdz = fdz/dz
!
  return
  end
!
!    =================================================
	subroutine DYCOMS_refine (nnz, dzz, fdz)
!    =================================================
!
USE gridno
USE shared_data
!
IMPLICIT NONE
!
      integer  :: k, ki, kf, nzmin, nnz
      real, dimension(1:nnz) :: fdz
      real, dimension(1:nnz+1) :: zw
      real  :: dzz
!
      real  :: zi, zf, zibar, zfbar, dzmin, fzi, domain_H
      real  :: b, d, dd, ddz, u, fmin, zz, zbar
      real  :: const_PI
      
      zi = 800.
      zf = 925.
      domain_H = 1500.
      dzmin = 5.0 
      fzi = 0.55
      
      ddz = domain_H / nnz

!     ... possibly adjust dzmin to fit evenly between zi and zf
      if ( zi < zf ) then
        nzmin = nint((zf - zi)/dzmin)
        dzmin = (zf - zi)/nzmin
      endif

!     ... indices of zi and zf
      ki = nint(fzi*nnz)
      kf = ki + (zf - zi)/dzmin

!     ... normalized zi and zf 
      zibar = ddz*ki
      zfbar = ddz*kf

!     ... ratio of minimum dz to H/nz
      fmin = dzmin / ddz

!     ... derived parameters
      const_PI = 2*acos(0.)
      b = (2.0)*( (zi/zibar) - fmin )
      d = (3.0/((domain_H - zfbar)**2) )*( (domain_H - zf)/(domain_H - zfbar) - fmin )

!     ... impossible grids
      if ( kf < ki .or. kf > nnz+1 ) then
         print *, 'Error :: ki, kf, nz = ', ki, kf, nnz
         stop
      endif
      if ( b < 0 ) then
         print *, 'Error :: b less than zero, b = ', b
         stop
      endif
      if ( d < 0 ) then
         print *, 'Error :: d less than zero, d = ', d
         stop
      endif

!     ... vertical positions of cell boundaries
      do k = 1, nnz+1
         zbar = ddz*(k-1)
         if ( k <= ki+1 ) then
            u  = const_PI*zbar/zibar
            zz = fmin*zbar+ (b*zibar/const_PI)*(u/2 - sin( 2.0*u )/4.0 )
         elseif ( k <= kf+1 ) then
            zz = zz + dzmin
         else
            dd = zbar - zfbar
            zz = zf + fmin*dd + d*dd*dd*dd/3.0
         endif
         zw(k) = zz
      enddo
      
!  Reevaluate fdz0
      do k = 1, nnz
        fdz(k) = (zw(k+1) - zw(k))/ddz
      enddo
!
      dzz = ddz
!
return
end

!    =================================================
 subroutine DYCOMS_refine_b (nnz, dzz, fdz)
!    =================================================
!
USE gridno
USE shared_data
!
IMPLICIT NONE
!
    integer  :: k, ki, kf, nzmin, nnz
    real, dimension(1:nnz) :: fdz
    real, dimension(1:nnz+1) :: zw
    real  :: dzz
!
    real  :: zi, zf, zibar, zfbar, dzmin, fzi, domain_H
    real  :: b, d, dd, ddz, u, fmin, zz, zbar
    real  :: const_PI

    zi = 750.
    zf = 925.
    domain_H = 3000.
    dzmin = 4.1 
    fzi = 0.4
    
    ddz = domain_H / nnz

!     ... possibly adjust dzmin to fit evenly between zi and zf
    if ( zi < zf ) then
      nzmin = nint((zf - zi)/dzmin)
      dzmin = (zf - zi)/nzmin
    endif

!     ... indices of zi and zf
    ki = nint(fzi*nnz)
    kf = ki + (zf - zi)/dzmin

!     ... normalized zi and zf 
    zibar = ddz*ki
    zfbar = ddz*kf

!     ... ratio of minimum dz to H/nz
    fmin = dzmin / ddz

!     ... derived parameters
    const_PI = 2*acos(0.)
    b = (2.0)*( (zi/zibar) - fmin )
    d = (3.0/((domain_H - zfbar)**2) )*( (domain_H - zf)/(domain_H - zfbar) - fmin )

!     ... impossible grids
    if ( kf < ki .or. kf > nnz+1 ) then
       print *, 'Error :: ki, kf, nz = ', ki, kf, nnz
       stop
    endif
    if ( b < 0 ) then
       print *, 'Error :: b less than zero, b = ', b
       stop
    endif
    if ( d < 0 ) then
       print *, 'Error :: d less than zero, d = ', d
       stop
    endif

!     ... vertical positions of cell boundaries
    do k = 1, nnz+1
       zbar = ddz*(k-1)
       if ( k <= ki+1 ) then
          u  = const_PI*zbar/zibar
          zz = fmin*zbar+ (b*zibar/const_PI)*(u/2 - sin( 2.0*u )/4.0 )
       elseif ( k <= kf+1 ) then
          zz = zz + dzmin
       else
          dd = zbar - zfbar
          zz = zf + fmin*dd + d*dd*dd*dd/3.0
       endif
       zw(k) = zz
    enddo
    
!  Reevaluate fdz0
    do k = 1, nnz
      fdz(k) = (zw(k+1) - zw(k))/ddz
    enddo
!
    dzz = ddz
!
return
end

!   =================================================
subroutine SE_Atlantic_refine (nnz, ddz, fdz)
  !   =================================================
  !
  use gridno
  use shared_data
    
  IMPLICIT NONE
  integer:: k,nnz
  !integer,parameter::  nnz=260
  real ::ddz,zz,h_res
  real :: domain_H 
  real :: domain_hr_H             !top of domain's part with higher resolution
  real,dimension(1:nnz):: fdz
  real,dimension(1:nnz+1):: zw
  
  domain_H = 6000.
  domain_hr_H = 2500.
  h_res=10
  ddz= domain_H/(nnz-1)
  
  k=1
  zz=0
  zw(k)=zz
  
  do while (zz<domain_hr_H)
    zz = zz + h_res
    k=k+1
    zw(k)= zz
  end do
  
  !do while (zz<=domain_H .and. k<=nnz+1)
  do while (k<nnz+1)	
    zz = zz+0.1*zz
    k=k+1
    zw(k) =zz
    !k=k+1
    !print*,k
  end do
  
  !  Reevaluate fdz0
  do k = 1, nnz
    fdz(k) = (zw(k+1) - zw(k))/ddz
  enddo
  
  return
end
