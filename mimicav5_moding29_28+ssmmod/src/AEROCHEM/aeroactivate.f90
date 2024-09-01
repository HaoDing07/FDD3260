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

MODULE aeroactivate

USE shared_aerosol_new

IMPLICIT NONE

PRIVATE

PUBLIC  :: activate

CONTAINS

subroutine activate ( s, s1, rg, sigma, den, A, kappa, nact, mact, rwet, ract )
  !
  REAL, INTENT(IN)  :: s                      !< actual supersatuation
  REAL, INTENT(IN)  :: s1                     !< supersatuation of previous substep
  REAL, INTENT(IN)  :: rg                     !< geometic mean radius
  REAL, INTENT(IN)  :: sigma                  !< standarddeviation aerosol number distribution
  REAL, INTENT(IN)  :: den                    !< density of aerosol population
  REAL, INTENT(IN)  :: A                      !< Kelvin curvature parameter
  REAL, INTENT(IN)  :: kappa                  !< Hygroskopicity aerosols
  !
  REAL, INTENT(OUT) :: nact                   !< activated number fraction
  REAL, INTENT(OUT) :: mact                   !< mean mass of activated fraction
  REAL, INTENT(OUT), OPTIONAL :: rwet         !< radius of wet activated fraction
  REAL, INTENT(OUT), OPTIONAL :: ract         !< radius of dry activated fraction
  !
  REAL  :: mg                                 !< geometic mean mass
  REAL  :: r                                  !< critical radius
  REAL  :: sigmam                             !< standarddeviation aerosol mass distribution
  REAL  :: rm                                 !< critical mass
  REAL  :: sigmar                             !< standarddeviation aerosol wet radius distribution
  REAL  :: rwg                                !< geometic mean wet radius 
  REAL  :: rr                                 !< critical wet radius
  REAL  :: r2                                 !< critical radius previous substep
  REAL  :: rr2                                !< critical wet radius previous substep
  REAL  :: x1, x2, y1, y2                     !< arguments of error function
  REAL  :: dn                                 !< Fraction of activated particles
!
!  Calculate activated fraction and mean mass of activated particles
!
  ! Calculate critical radius and activated fraction
  r = A/3. * (4./kappa/s**2.)**(1./3.)
  x1 = log(r/rg) / (sqrt(2.)*log(sigma))

  ! Critical radius of previously activated particles
  r2 = A/3. * (4./kappa/max(s1,1.e-10)**2.)**(1./3.)
  x2 = log(r2/rg) / (sqrt(2.)*log(sigma))

  ! Fraction of particles activated in the current iteration
  dn = 0.5*(erf(x2) - erf(x1))
!
  if ( dn > 1.e-8 ) then
    ! Parameters of aerosol mass distribution (by definition, also log-normal)
    rm = 4./3.*pi*den*r**3.
    mg = 4./3.*pi*den*rg**3.
    sigmam = sigma**3.
!
    ! Activated fraction
    nact = 0.5*(1. - erf(x1))
!
    ! Calculate mean mass of activated particles
    mact = mg*exp(0.5*log(sigmam)*log(sigmam)) *   &
        0.5*(1. - erf((log(rm/mg) - log(sigmam)*log(sigmam)) / (sqrt(2.)*log(sigmam))))
!
!  If required, calculate mean wet and dry radii of activated particles
!  What we want here is the size of particles activated between two iterations
!
    ract = 0.
    rwet = 1.e9
    if ( present(ract) .and. present(rwet) ) then
      ! Parameters of wet size distribution (by definition, also log-normal)
      rwg = sqrt(3.*kappa*rg**3./A)
      sigmar = sigma**1.5
      rr = sqrt(3.*kappa*r**3./A)
      rr2 = sqrt(3.*kappa*r2**3./A)
!
      ! Calculate mean dry radius of activated particles
      y1 = (log(r/rg) - log(sigma)*log(sigma)) / (sqrt(2.)*log(sigma))
      y2 = (log(r2/rg) - log(sigma)*log(sigma)) / (sqrt(2.)*log(sigma))
      ract = rg*exp(0.5*log(sigma)*log(sigma)) * (erf(y1) - erf(y2)) / (erf(x1) - erf(x2))
!
      ! Calculate mean wet radius of activated particles
      y1 = (log(rr/rwg) - log(sigmar)*log(sigmar)) / (sqrt(2.)*log(sigmar))
      y2 = (log(rr2/rwg) - log(sigmar)*log(sigmar)) / (sqrt(2.)*log(sigmar))
      rwet = rwg*exp(0.5*log(sigmar)*log(sigmar)) * (erf(y1) - erf(y2)) / (erf(x1) - erf(x2))
    endif
  else
    nact = 0.
    mact = 0.
    ract = 0.
    rwet = 1.e9
  endif
!  
  return
END subroutine activate

END MODULE aeroactivate

