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
!  SHARED_HYDRO.F:                   
!
!  Purpose:
!	Module containing aerosol variables (replaces aerodef.h, aeropar.h) 
!
!  Author
!	Julien Savre, MISU
!
! ================================================================

	module shared_aerosol_new

	use gridno
	use shared_data
	use typedef_aerosol_new

	SAVE

        type(aero_flg_list) :: aero_flg
        real, pointer, dimension(:,:,:)  :: tab_coag
	type (aeroprop), pointer, dimension(:) :: aero
	type (aero_1d), pointer, dimension(:)  :: aero0
	type (aero_3d), pointer, dimension(:)  :: aero3d2, aero3d, aero3ds, aero3d_mean

        type(aero_init), dimension(30) :: aeroi
	type(aero_init) :: aeror

	!
	! --- aerosol fixed parameters, constants...
	!
	! 1) sulfate, 2) BC, 3) sea salt, 4) orga1
	type(elementary), dimension(1:nelem) 				   &
		 :: elem = (/						   &
! 					   (kappa, rho, mw, name),
		 		elementary(0.55, 1840., 98., 'SULF'),      &
				elementary(0.20, 600., 12., 'BC'),	   &                   ! Before kappa= 0.01
		 		elementary(1.12, 2180., 58.4, 'SALT'),     &               
				elementary(0.20, 1560., 104., 'ORGA') /)                   ! Before kappa = 0.06
!
	type(aero_param)		  						   &
		 :: aerop = 								   &
!				 (critn, nucint, aitint, bcagemass,						    
				aero_param(100., 3.8e-3, 23.8e-3, 8.e5)	
!
	type(aero_bulk)		  							   &
		 :: aerob(2) = (/ 							   &
!				 (mmol, rho, nuphi, beta, b,	    
				aero_bulk(0.0585, 2160., 2., 0.5, 1.33),		   &
				aero_bulk(0.132, 1770., 2.3, 0.5, 0.55)	  	 /)
!
!  List of elementary functions
!
	public :: initialize_mix
	public :: have_element
	public :: cal_rgeo
	public :: cal_rav
	public :: cal_mav
	public :: cal_lnmom
	public :: cal_logmom
	public :: cal_inclogn
        public :: calerf
        public :: setup_aero_flags

        CONTAINS
  
! ============================
	
  subroutine initialize_mix (i)
  
  integer :: i, j, n

#ifdef AERO_ENABLE    
  aero(i)%mix%rho   = 0.0
  aero(i)%mix%kappa = 0.0
  aero(i)%mix%mw    = 0.0

  n = 0
  do j = 1, 4
    if (aero(i)%init%present(j)) then
      n = n + 1
      aero(i)%mix%kappa = aero(i)%mix%kappa + aero(i)%init%frac(j) * elem(j)%kappa/elem(j)%rho
      aero(i)%mix%mw = aero(i)%mix%mw + aero(i)%init%frac(j) * elem(j)%mw
      aero(i)%mix%rho = aero(i)%mix%rho + aero(i)%init%frac(j) * elem(j)%rho
    endif
  enddo
  aero(i)%mix%kappa = aero(i)%mix%kappa*aero(i)%mix%rho
#endif

  end subroutine initialize_mix
  
! ============================
	
  subroutine have_element (im, name, present, ielem)
  
  logical :: present
  character(len=4) :: name
  integer :: im, ielem

#ifdef AERO_ENABLE
  integer :: iel
  
  if (name == 'SULF') then
    iel = 1
  else if (name == 'BC  ') then
    iel = 2  
  else if (name == 'SALT') then
    iel = 3  
  else if (name == 'ORGA') then
    iel = 4  
  endif
  
  present = .false.
  ielem = 0
  if ( aero(im)%init%present(iel) ) then
    present=.true.
    ielem = iel
  endif
#endif
  
  end subroutine have_element
	
! ============================

  real function cal_rgeo (s, rhop, mm, nn)
  
  real :: s, mm, nn, rhop, raa
!
  if (nn > 1.e-3) then
    raa = (mm/nn / (4.*pi/3.*rhop))**(1./3.)
  else
    raa = 0.
  endif
!
!  Calculate modal diameter
!  
  if (s /= 0.) then
    cal_rgeo = raa * exp(-1.5*log(s)*log(s))
  else
    cal_rgeo = raa
  endif  
  
  return
  end function
	
! ============================

  real function cal_rav (s, rhop, mm, nn)
  
  real :: s, mm, nn, rhop, rr
!
  rr = cal_rgeo(s, rhop, mm, nn)
!
  cal_rav = rr * exp(0.5*log(s)*log(s))
!
  return
  end function

	
! ============================

  real function cal_mav (s, rhop, mm, nn)
  
  real :: s, mm, nn, rhop, rr
!
  rr = cal_rgeo(s, rhop, mm, nn)
!
  cal_mav = 4./3.*pi*rhop*rr**3. * exp(4.5*log(s)*log(s))
!
  return
  end function

! ============================

  real function cal_lnmom (n, s, d)
  
  real :: n, s, d, mu
!
!  Calculate nth moment, with d the modal diameter
!
  if (s /= 0.) then
    cal_lnmom = d**n * exp( 0.5*(n*log(s))**2. )
  else
    cal_lnmom = d**n
  endif

  return
  end function

! ============================

  function cal_logmom (i,s,m,n,rho)

  real     :: i, m, n, s, rho
  real     :: d, cal_logmom

!  Calculate ith moment, from mass and number

  d = 2.*cal_rgeo (s, rho, m, n)
!  
  cal_logmom = cal_lnmom (i, s, d)

  return
  end function cal_logmom

! ================================================

  real function cal_inclogn (mu,sigma,phi)

  real     :: mu, sigma, phi
  real     :: xx, s, erf1

  xx = 2.*phi - 1.
  s  = sign(1.,xx)
  erf1 = cal_inverror (min(abs(xx),1.-1.e-8))

  xx = log(mu) + sqrt(2.)*log(sigma)*erf1*s
  cal_inclogn = exp(xx)

  return
  end function
      
! ================================================

  subroutine setup_aero_flags ()

    if (laero < 0) then
      aero_flg%tran     = .false.
      aero_flg%act      = .false.
      aero_flg%act_scv  = .false.
      aero_flg%reg      = .false.
      aero_flg%imp_scv  = .false.
      aero_flg%chem     = .false.
    elseif (laero .eq. 0) then
      aero_flg%tran     = .false.
      aero_flg%act      = .true.
      aero_flg%act_scv  = .false.
      aero_flg%reg      = .false.
      aero_flg%imp_scv  = .false.
      aero_flg%chem     = .false.
    elseif  (laero .eq. 1) then
      aero_flg%tran     = .true.
      aero_flg%act      = .true.
      aero_flg%act_scv  = .true.
      aero_flg%reg      = .false.
      aero_flg%imp_scv  = .false.
      aero_flg%chem     = .false.
    elseif  (laero .eq. 2) then
      aero_flg%tran     = .true.
      aero_flg%act      = .true.
      aero_flg%act_scv  = .true.
      aero_flg%reg      = .true.
      aero_flg%imp_scv  = .false.
      aero_flg%chem     = .false.
    elseif  (laero .eq. 3) then
      aero_flg%tran     = .true.
      aero_flg%act      = .true.
      aero_flg%act_scv  = .true.
      aero_flg%reg      = .true.
      aero_flg%imp_scv  = .true.
      aero_flg%chem     = .false.
    elseif  (laero > 3) then
      !aero_flg set directly in cm.nml
    endif

    aero_flg%any = aero_flg%tran .or. aero_flg%act .or. aero_flg%act_scv .or. aero_flg%imp_scv .or. aero_flg%reg .or. aero_flg%chem

  end subroutine setup_aero_flags
      		
! ============================

  function calerf ( arg, jint )

!*****************************************************************************80
!
!! CALERF computes various forms of the error function.
!
!  Discussion:
!
!    This routine evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
!    for a real argument x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    Rational Chebyshev Approximations for the Error Function,
!    Mathematics of Computation,
!    Volume 23, Number 107, July 1969, pages 631-638.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT is 1, the
!    argument must be less than XBIG.  If JINT is 2, the argument
!    must lie between XNEG and XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = erf(x);
!    1, RESULT = erfc(x) = 1 - erf(x);
!    2, RESULT = exp(x*x)*erfc(x) = exp(x*x) - erf(x*x)*erf(x).
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    0, erf(x);
!    1, erfc(x);
!    2, exp(x*x)*erfc(x).
!
  implicit none

  real a(5)
  real arg
  real b(4)
  real c(9)
  real d(8)
  real del
  integer i
  integer jint
  real p(6)
  real q(5)
  real result
  real sixten
  real sqrpi
  real thresh
  real x
  real xbig
  real xden
  real xhuge
  real xinf
  real xmax
  real xneg
  real xnum
  real xsmall
  real y
  real ysq
  real calerf
!
!  Mathematical constants
!
  data sqrpi / 5.6418958354775628695d-1 /
  data thresh / 0.46875d0 /
  data sixten / 16.0d0 /
!
!  Machine-dependent constants
!
  data xneg / -26.628d0 /
  data xsmall /1.11d-16/
  data xbig /26.543d0 /
  data xhuge /6.71d7/
  
  xmax=huge(x)
  xinf=sqrt(huge(x))
!
!  Coefficients for approximation to  erf  in first interval
!
  data a/3.16112374387056560d00,1.13864154151050156d02, &
         3.77485237685302021d02,3.20937758913846947d03, &
         1.85777706184603153d-1/
  data b/2.36012909523441209d01,2.44024637934444173d02, &
         1.28261652607737228d03,2.84423683343917062d03/
!
!  Coefficients for approximation to  erfc  in second interval
!
  data c/5.64188496988670089d-1,8.88314979438837594d0, &
         6.61191906371416295d01,2.98635138197400131d02, &
         8.81952221241769090d02,1.71204761263407058d03, &
         2.05107837782607147d03,1.23033935479799725d03, &
         2.15311535474403846d-8/
  data d/1.57449261107098347d01,1.17693950891312499d02, &
         5.37181101862009858d02,1.62138957456669019d03, &
         3.29079923573345963d03,4.36261909014324716d03, &
         3.43936767414372164d03,1.23033935480374942d03/
!
!  Coefficients for approximation to  erfc  in third interval
!
  data p/3.05326634961232344d-1,3.60344899949804439d-1, &
         1.25781726111229246d-1,1.60837851487422766d-2, &
         6.58749161529837803d-4,1.63153871373020978d-2/
  data q/2.56852019228982242d00,1.87295284992346047d00, &
         5.27905102951428412d-1,6.05183413124413191d-2, &
         2.33520497626869185d-3/

  x = arg
  y = abs ( x )
!
!  Evaluate erf for |X| <= 0.46875.
!
  if ( y <= thresh ) then

    result = -1.0D+00

    ysq = 0.0D+00
    if ( xsmall < y ) then
      ysq = y * y
    end if

    xnum = a(5) * ysq
    xden = ysq

    do i = 1, 3
      xnum = ( xnum + a(i) ) * ysq
      xden = ( xden + b(i) ) * ysq
    end do

    result = x * ( xnum + a(4) ) / ( xden + b(4) )

    if ( jint /= 0 ) then
      result = 1.0D+00 - result
    end if

    if ( jint == 2 ) then
      result = exp ( ysq ) * result
    end if

    calerf = result

    return
!
!  Evaluate erfc for 0.46875 <= |X| <= 4.0.
!
   else if ( y <= 4.0D+00 ) then

     xnum = c(9) * y
     xden = y

     do i = 1, 7
       xnum = ( xnum + c(i) ) * y
       xden = ( xden + d(i) ) * y
     end do

     result = ( xnum + c(8) ) / ( xden + d(8) )

     if ( jint /= 2 ) then
       ysq = aint ( y * sixten ) / sixten
       del = ( y - ysq ) * ( y + ysq )
       result = exp ( -ysq * ysq ) * exp ( -del ) * result
     end if
!
!  Evaluate erfc for 4.0 < |X|.
!
   else

     result = 0.0D+00

     if ( xbig <= y ) then

       if ( jint /= 2 .or. xmax <= y ) then
         go to 300
       end if

       if ( xhuge <= y ) then
         result = sqrpi / y
         go to 300
       end if

     end if

     ysq = 1.0D+00 / ( y * y )
     xnum = p(6) * ysq
     xden = ysq
     do i = 1, 4
       xnum = ( xnum + p(i) ) * ysq
       xden = ( xden + q(i) ) * ysq
      end do

      result = ysq * ( xnum + p(5) ) / ( xden + q(5) )
      result = ( sqrpi -  result ) / y

      if ( jint /= 2 ) then
        ysq = aint ( y * sixten ) / sixten
        del = ( y - ysq ) * ( y + ysq )
        result = exp ( -ysq * ysq ) * exp ( -del ) * result
      end if

  end if
!
!  Fix up for negative argument, erf, etc.
!
  300 continue

  if ( jint == 0 ) then

    result = ( 0.5D+00 - result ) + 0.5D+00
    if ( x < 0.0D+00 ) then
      result = -result
    end if

  else if ( jint == 1 ) then

    if ( x < 0.0D+00 ) then
      result = 2.0D+00 - result
    end if

  else

    if ( x < 0.0D+00 ) then

      if ( x < xneg ) then
        result = xinf
      else
        ysq = aint ( x * sixten ) / sixten
        del = ( x - ysq ) * ( x + ysq )
        y = exp ( ysq * ysq ) * exp ( del )
        result = ( y + y ) - result
      end if

    end if

  end if
  
  calerf = result

  return
  end function

! ============================
	
      real function cal_inverror (p)
      
      IMPLICIT NONE
      
      real :: p,p_low,p_high
      real :: a1,a2,a3,a4,a5,a6
      real :: b1,b2,b3,b4,b5
      real :: c1,c2,c3,c4,c5,c6
      real :: d1,d2,d3,d4
      real :: z,q,r
      
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1.-p_low
      
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
      
201   q=sqrt(-2.*log(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/			&
     	((((d1*q+d2)*q+d3)*q+d4)*q+1.)
      goto 204
      
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
      
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/		&
     	(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
      
302   if((p.gt.p_high).and.(p.lt.1.)) goto 203
203   q=sqrt(-2.*log(1.-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/			&
     	 ((((d1*q+d2)*q+d3)*q+d4)*q+1.)

204   cal_inverror=z
      
      return
      end function				

  end module shared_aerosol_new
