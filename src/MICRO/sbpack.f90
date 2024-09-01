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
!  MICROPACK1:
!	Package of subroutines related to microphysics:
!
!       Warm phase microphysics: treated by default using Seifert & Beheng's scheme
!       aucr (droplet autoconversion and droplet self-collection), 
!       ccr (cloud drops accretion), scr (rain self-collection & breakup)
!
!  Author:
!	Julien Savre, MISU
!
! ================================================================
!
  module sbpack
!
  USE shared_data
  USE shared_hydro
  USE shared_thermo
  USE shared_pressure
!
  IMPLICIT NONE
!
  private
!
  real, parameter :: drim = 500.e-6, qrim = 0.1e-3, deq = 0.9e-3, dbu = 0.35e-3, dmin = 10.e-6, dhm = 50.e-6
  real, parameter :: nus = 1., kc = 9.44e9, kr = 5.25, krr = 7.12, kapr = 60.7, hmr = 3.5e8, kliu = 1.08e10
  real, parameter :: alphai = 0.68, alphas = 0.01
  real, parameter :: ea0 = 0.8
  real, parameter :: krac = 1.e-3
  real, parameter :: Re_max = 10000.
!
  public :: aucr, ccr, scr, scc, kess, liu_daum
  public :: sbevap, sbcond, sbdep, sbmelt, sbcollec, sbrime, sbrimemelt, hallett_mossop
  public :: fcaltdi, cshed
!
  contains
!
!  ==================================================
   subroutine sbcond (k, h, qs, es, Tem, xfkd, xfmu, q, n, cd)
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Condensation of water drops, with ventilation
!
! ================================================================
!
  integer  :: k, h
  real     :: q, n, d, m
  real     :: Tem, es, qs
  real     :: xfkd, xfmu, Fv, ff
  real     :: cd
!
!----------------------------------------------------------!
! 
  ff = 1. / (1./(den0(k)*xfkd) + qs*cal_flv(Tem)/(Kt*Tem)*(cal_flv(Tem)/(Rw*Tem) - 1.))
!
  m = cal_avm(h, q, n)
  d = hydrop(h)%cm*m**(hydrop(h)%am) 
!
  Fv = cal_fvq (h,den0(k),xfmu,xfkd,q,n,pxx(k,h)) 
!
  cd = 2.*pi * ff * Fv * d * n
  cd = max(cd,0.0)
!
return
end
!
!  ==================================================
   subroutine sbdep (k, h, qsi, esi, Tem, xfkd, xfmu, q, n, dep)
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Deposition of ice, with ventilation
!
! ================================================================
!
IMPLICIT NONE
!
  integer  :: k, h
  real     :: n, q, d, m
  real     :: Tem, esi, qsi
  real     :: dep
  real     :: xfkd, xfmu, ff, Fv
!
!----------------------------------------------------------!
!
  ff = 1. / (1./(den0(k)*xfkd) + qsi*cal_fls(Tem)/(Kt*Tem)*(cal_fls(Tem)/(Rw*Tem) - 1.))
!
  m = cal_avm(h, q, n)
  d = hydrop(h)%cm*m**(hydrop(h)%am)
!
  Fv = cal_fvq (h,den0(k),xfmu,xfkd,q,n,pxx(k,h))
!
  dep = 4.*pi*hydrop(h)%cap * ff * Fv * d * n
  dep = max(dep,0.0)
!
return
end
!
!  ==================================================
   subroutine sbmelt (k,h,q,n,Tem,esw,es0,xfkd,xfmu,xmr,xmrn,xmin)
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Melting of ice to rain
!
! ================================================================
!
  integer  :: k, h
  real     :: q, n, m, Tem
  real     :: xfkd, xfmu, esw, es0
  real     :: xmr, xmrn, xmin
  real     :: ff, Fv, a, d, m1
!
!----------------------------------------------------------!  
!
  !ff = (Kt*Kt/(den0(k)*cp_a*xfkd)*(Tem - 273.15) + cal_fls(Tem)*xfkd/Rw*(esw/Tem - es0/273.15)) / cal_flm(Tem)
  ff = (Kt*Kt/(den0(k)*cp_a*xfkd)*(Tem - 273.15)) / cal_flm(Tem)
!
  m = cal_avm(h, q, n)
  d = hydrop(h)%cm*m**(hydrop(h)%am) 
!
  Fv = cal_fvq (h,den0(k),xfmu,xfkd,q,n,pxx(k,h)) 
!
  xmr = 2.*pi * ff * Fv * n * d
  xmr = max(min(xmr,q/dt0),0.)
!  
  if ( moments == 2 ) then
    a = cal_fvn (h,den0(k),xfmu,xfkd,q,n,pxx(k,h)) / cal_fvq (h,den0(k),xfmu,xfkd,q,n,pxx(k,h)) 
    xmin = max(min(a*xmr/m,n/dt0),0.)
!
    if ( h /= hail ) then
      xmrn = xmin
    else
      m1 = (1.e-3/hydrop(h)%cm)**(1./hydrop(h)%am)
      xmrn = xmr/m1
    endif
  endif
!
return
end
!
!  ==================================================
   subroutine sbevap (k, h, q, n, xfkd, xfmu, dq, dn)
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Evaporation of cloud drops
!
! ================================================================
!
  integer  :: h, k
  real     :: xfkd, xfmu
  real     :: q, n, m, a, dq, dn
!
!----------------------------------------------------------!
!  
  if ( h > drop .and. dq < 0. .and. q > 1.e-10 ) then
    a = cal_fvn (h,den0(k),xfmu,xfkd,q,n,pxx(k,h)) / cal_fvq (h,den0(k),xfmu,xfkd,q,n,pxx(k,h)) 
    m = cal_avm(h, q, n)
    dn = max(a*dq/m,-n/dt0)
  else if ( (lndrop == 1 .and. dq < 0. .and. q < 1.e-10) .or. (h > drop .and. dq < 0. .and. q < 1.e-10) ) then
    dn = -n/dt0
  endif
!
return
end
!
!  ==================================================
   subroutine aucr (k,dens,qc,qr,nc,nr,auq,aunc,aunr)
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Autoconversion from cloud to rain + droplets self-collection
!
! ================================================================
!
  integer  :: k
  real     :: qc, qr, nc, nr
  real     :: dens
  real     :: auq, aunc,aunr
  real     :: xx, tau, fautau, qcdens, qrdens, ncdens, nrdens
!
!  Conversions
!
  xx = cal_avm(drop, qc, nc)
  qcdens=qc*dens
  qrdens=qr*dens
  ncdens=nc*dens 
  nrdens=nr*dens 
!
!  Calculate model parameters
!
  if ( ncdens > xnc_min .and. qcdens > qc_min ) then
    tau = max(min(1. - qcdens/(qcdens + qrdens),0.99),0.)
    fautau = max(400.*tau**0.7*(1 - tau**0.7)**3., 0.)
!
!  Variation of rain mixing ratio through autoconversion
!
    auq = kc/(20.*xauto) * (nus + 4.)*(nus + 2.)/(nus + 1.)**2. * qcdens**2. * xx**2. * (1. + fautau/(1. - tau)**2.) * pxx(k,drop)**2.
    auq = min(max(auq/dens,0.0),qc/dt0)
!  
! Variation of cloud droplet number: self-collection + auto-conversion
!
    if (moments == 2) then
      aunc = min(max(0.0,auq/xx),nc/dt0)
      aunr = min(max(0.0,auq/xauto),(nc-nr)/dt0)
    endif
  endif
!
return
end
!
!  ==================================================  
  subroutine kess ( dens, qc, nc, auq, aunc, aunr )
!  ==================================================

    real, intent(in) :: dens, qc, nc
    real, intent(out) :: auq, aunc, aunr
!
!----------------------------------------------------------!
!
    if (qc*dens > qc_min) then
      auq = max(krac*(qc - qauto), 0.)
      aunc = nc/qc*auq
      aunr = auq/xauto
    endif
! 
  end subroutine kess
!
!  ==================================================  
  subroutine liu_daum ( k, dens, qc, nc, nr, auq, aunc, aunr )
!  ==================================================

    integer, intent(in) :: k
    real, intent(in) :: dens, qc, nc, nr
    real, intent(out) :: auq, aunc, aunr
    real :: qcdens, rauto, xx, r6, beta6
!
!----------------------------------------------------------!
!
    xx = cal_avm(drop, qc, nc)
    qcdens = qc*dens
!
    if ( qc*dens > qc_min .and. nc*dens > xnc_min ) then
      beta6 = hydroc(drop)%cr6
      rauto = beta6**(1./6.) * 0.5*hydrop(drop)%cm*(xauto)**hydrop(drop)%am
      r6 = beta6**(1./6.) * 0.5*hydrop(drop)%cm*(xx)**hydrop(drop)%am
!
      if ( r6 >= rauto ) then
        auq = kliu * beta6 * qcdens**2. * xx * pxx(k,drop)**2.
        auq = min(max(auq/dens,0.0),qc/dt0)
!
        if (moments == 2) then
          aunc = min(max(0.0,auq/xx),nc/dt0)
          aunr = min(max(0.0,auq/xauto),(nc-nr)/dt0)
        endif
      endif
    endif
! 
  end subroutine liu_daum
!
!  ==================================================
   subroutine ccr (k,dens,qc,qr,nc,nr,clr,clrn)
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Collection of cloud droplets by rain (accretion)
!
! ================================================================
!
  integer  :: k
  real     :: qc, qr, nc, nr, m
  real     :: clrn, clr
  real     :: tau, factau, qcdens, qrdens, dens
!
!----------------------------------------------------------!
!                     Seifert model                        !
!----------------------------------------------------------!
!
!  Conversions
!
          qcdens=qc*dens
          qrdens=qr*dens
!
!  Calculate model parameters
!
  	  if (qrdens > qr_min .and. qcdens > qc_min) then
            tau = max(min(1. - qcdens/(qcdens + qrdens),1.),0.)
            factau = (tau / (tau + 5.E-4))**4.
!
!  Variation of rain/drop mixing ratio through accretion
!
            clr = kr * qcdens * qrdens * factau * pxx(k,rain)
            clr = min(max(clr/dens,0.0),qc/dt0)
!
!  Variation of droplet number through accretion
!
	    if (moments == 2 .and. lndrop == 1) then
              clrn = min(max(nc/qc*clr,0.0),nc/dt0)
	    endif
          endif 
!
return
end
!
!  ==================================================
   subroutine scr (k,dens,qr,nr,crrn)
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Self collection of rain & breakup
!
! ================================================================
!
  integer  :: k
  real     :: nr, qr, qrdens, nrdens, lambdar, crrn
  real     :: dens, phibr, dm
!
!----------------------------------------------------------!
!                     Seifert model                        !
!----------------------------------------------------------!
!
!  Conversions
!
  qrdens=qr*dens
  nrdens=nr*dens 
  lambdar = cal_lambda(rain, qr, nr)
!
!  Add self collection of rain droplets
!
  crrn = krr * nrdens * qrdens * (1. + kapr/lambdar)**(-9.) * pxx(k,rain)
  crrn = max(crrn/dens,0.0)
!
!  Breakup function from Seifert & Beheng
! 
  dm  = cal_avd(rain,qr,nr)
  if (dm < dbu) then
    phibr = 1.
  else if (dm >= dbu) then
    phibr = -900.*(dm - deq)
  endif
!
!  Total self-collection of rain:
!
  crrn = phibr*crrn
!
return
end
!
!  ==================================================
   subroutine scc (k,dens,qc,nc,cccn)
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Self collection of cloud drops
!
! ================================================================
!
  integer  :: k
  real     :: nc, qc, qcdens, ncdens, dens
  real     :: cccn
!
!----------------------------------------------------------!
!                     Seifert model                        !
!----------------------------------------------------------!
!
!  Conversions
!
  qcdens=qc*dens
  ncdens=nc*dens 
!
!  Self collection of cloud droplets (no breakup)
!
  cccn = kc * (nus+2) / (nus+1) * qcdens**2. * pxx(k,drop)**2.
  cccn = min(max(cccn/dens,0.0),nc/dt0)
!
return
end
!
!  ==================================================
   subroutine sbcollec (i,j,dens,Tem,qi,ni,qj,nj,cli,clj,cln)
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Generalized collection process involving ice particles
!	Following Seifert & Beheng
!
! ================================================================
!
  integer  :: i,j
  real     :: dens, Tem, qi, qj, ni, nj, di, dj, vi, vj, mi, mj
  real     :: qidens, qjdens, nidens, njdens, dvqij, dvqji, dvn, intqij, intqji, intn, ea, es
  real     :: cli, clj, cln
!
!----------------------------------------------------------!  
!
!  Conversions
!
  qidens=qi*dens
  nidens=ni*dens
  qjdens=qj*dens
  njdens=nj*dens
!
  mi = cal_avm(i, qi, ni)
  mj = cal_avm(j, qj, nj)
!
  di = hydrop(i)%cm*mi**(hydrop(i)%am)
  dj = hydrop(j)%cm*mj**(hydrop(j)%am)
  vi = hydrop(i)%ctv*mi**(hydrop(i)%btv)
  vj = hydrop(j)%ctv*mj**(hydrop(j)%btv)
!
  if (j == drop) then
    if (cal_avd(j,qj,nj) < 20.e-6) then
      ea = 0.
    else if ( cal_avd(j,qj,nj) >= 20.e-6 .and. cal_avd(j,qj,nj) < 100.e-6 ) then
      ea = (cal_avd(j,qj,nj) - 20.e-6)/(100.e-6 - 20.e-6)
    else
      ea = 1.
    endif
    if ( cal_avd(i,qi,ni) < 200.e-6 ) then
      ea = 0.
    else
      if (i == snow .or. i == ice) ea = ea0*ea
    endif
  else if (i == rain .or. j == rain) then
    ea = 1.
  else
!    ea = 0.1	! Rutledge and Hobbs 
    ea = exp(0.09*(Tem - 273.15))	! Original SB scheme
  endif
!
!  Q
!
  intqij = hydroc(i)%deltan*di*di + 2.*hydroc(i)%deltaijq(j)*di*dj + hydroc(j)%deltaq*dj*dj
  intqji = hydroc(j)%deltan*dj*dj + 2.*hydroc(j)%deltaijq(i)*di*dj + hydroc(i)%deltaq*di*di
!
  dvqij = sqrt( hydroc(i)%thetan*vi*vi - 2.*hydroc(i)%thetaijq(j)*vi*vj + hydroc(j)%thetaq*vj*vj )
  dvqji = sqrt( hydroc(j)%thetan*vj*vj - 2.*hydroc(j)%thetaijq(i)*vi*vj + hydroc(i)%thetaq*vi*vi )
!
  cli = 0.25*pi*ea*qidens*njdens*intqji*dvqji
  clj = 0.25*pi*ea*qjdens*nidens*intqij*dvqij
!
  cli = min(max(0.0,cli/dens),qi/dt0)
  clj = min(max(0.0,clj/dens),qj/dt0)
!
!  N
!
  if ( moments == 2 ) then
    intn = hydroc(i)%deltan*di*di + 2.*hydroc(i)%deltaijn(j)*di*dj + hydroc(j)%deltan*dj*dj
    dvn = sqrt( hydroc(i)%thetan*vi*vi - 2.*hydroc(i)%thetaijn(j)*vi*vj + hydroc(j)%thetan*vj*vj )
!
    cln = 0.25*pi*ea*nidens*njdens*intn*dvn
    cln = min(max(0.0,cln/dens),nj/dt0)
  endif
!
return
end
!
!  ==================================================
   subroutine sbrime (h,Tem,qi,ni,frim)
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Generalized collection process involving ice particles
!	Following Seifert & Beheng
!
! ================================================================
!
  integer  :: h
  real     :: Tem, qi, ni
  real     :: alpha, di, qg, r_m, frim
!
!----------------------------------------------------------!  
!
  if (h == ice) then
    alpha = alphai*rhow/rhoi
  else if (h == snow) then
    alpha = alphas*rhow/rhoi
  endif
!
  di = cal_avd(h, qi, ni)
  qg = ni*pi/6.*rhoi*di**3.
!
!  Riming occurs only if T<0 and ice particles are large enough: follows Seifert & Beheng
!
  if (di >= drim .and. Tem < 273.15) then
    r_m = qg/qi - 1.
    frim = max(min(1./(alpha*r_m), 1.), 0.)
  else
    frim = 0.
  endif
!
  if (.not.lrime) frim = 0.
!
return
end
!
!  ==================================================
   subroutine sbrimemelt (Tem,riml,rimr,melt)
!  ==================================================
!
! ================================================================
!  Purpose:
!	Melt of graupel during riming
! ================================================================
!
  real :: Tem, rimr, riml, melt
!
!----------------------------------------------------------!  
!
  melt = cp_l/cal_flm(Tem)*(Tem - 273.15)*(rimr + riml)
!
return
end
!
!  ==================================================
   subroutine hallett_mossop ( Tem, qc, nc, rime, hm, hmn )
!  ==================================================
!
! ================================================================
!  Purpose:
! 	 Calculation of secondary ice multiplication, hallett-mossop process
! ================================================================
!
  real  :: Tem, qc, nc, dc, mc, rime, hm, hmn, ff, xmini
!
!----------------------------------------------------------!
!
  mc = cal_avm(drop, qc, nc)
  dc = hydrop(drop)%cm*mc**(hydrop(drop)%am)
  xmini = (dmin/hydrop(ice)%cm)**(1./hydrop(ice)%am)
!
  if (dc > dhm .and. Tem > 267.5 .and. Tem <= 270.) then
    ff = 1. - (Tem - 267.5) / (270. - 267.5)
  else if (dc > dhm .and. Tem >= 265. .and. Tem <= 267.5) then
    ff = (Tem - 265.) / (267.5 - 265.)
  else
    ff = 0.
  endif
!
  hm = hm + hmr*ff*rime*xmini
!
  if (moments == 2) hmn = hmn + hmr*ff*rime
!
!----------------------------------------------------------!
!
return
end
!
!  ==================================================
   subroutine fcaltdi ( k, h, T, ql, esi, xfkd, xfmu, q, n, dql, dqi, Td, cghr )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Calculation of surface temperature of particle h
!
! ================================================================
!
  integer :: k, h
  real    :: T, esi, xfkd, xfmu, ql, q, n, m, dql, dqi, Td, cghr
  real    :: ff, ffcrit, esc, Fv, d, I1, I2, I3, dTcrit, Tcrit
!
!----------------------------------------------------------!
!
  dTcrit = min(-0.5*(T - 273.15),3.)
  Tcrit = 273.15 - dTcrit
  esc = cal_esi( Tcrit )
!
  ffcrit = Kt*Kt/(den0(k)*cp_a*xfkd) + cal_fls(Tcrit)*xfkd*esi/(Rw*Tcrit*Tcrit)*(cal_fls(Tcrit)/(Rw*Tcrit) - 1.)
  ff = Kt*Kt/(den0(k)*cp_a*xfkd)*(T - Tcrit) + cal_fls(T)*xfkd/Rw*(esi/T - esc/Tcrit)
!
  Fv = cal_fvq (h,den0(k),xfmu,xfkd,q,n,pxx(k,h))
!
  m = cal_avm(h, q, n)
  d = hydrop(h)%cm*m**(hydrop(h)%am)
!
  I1 = cal_flm(T)*dql + 4.*pi*hydrop(h)%cap*ff*Fv*n*d + (cp_l*dql + cp_i*dqi)*(T - Tcrit)
  I2 = 4.*pi*hydrop(h)%cap*ffcrit*Fv*n*d + (cp_l*dql + cp_i*dqi)
!
  Td = I1 / I2
!
  if ( d > 1.e-3 .and. T > 248.15 .and. ql > max(-2.e-4*(T - 273.15),1.e-3) ) cghr = 0.368*exp(0.333*min(Td,dTcrit)) 	! or 1.
!
return
end
!
!  ==================================================
   subroutine cshed ( h, T, ql, q, n, dql, dqi, cghs, cghsn )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Shedding rate for wet graupel and hail particles
!		Shedding rate and threshold follow Garcia-Garcia and List, 1992
!		Shed liquid water is assumed to be in the form of 1mm rain drops
!
! ================================================================
!
  integer :: h
  real    :: T, q, n, dql, dqi, ql, qlc, m1
  real    :: cgh, cghs, cghsn
!
!----------------------------------------------------------!
!
  qlc = 1.e-3*(2.1 - 0.18*(T - 273.15))
!
  if ( ql > qlc ) then
    cgh = (dql + dqi) / (1. + (100. + 4.*(T - 273.15))*(ql - qlc))
  else
    cgh = dql + dqi 
  endif
!
  cghs = max(dql + dqi - cgh, 0.0)
!
  m1 = (1.e-3/hydrop(rain)%cm)**(1./hydrop(rain)%am)
  if (moments == 2) cghsn = cghs / m1
!
return
end
!
!  ==================================================
  function cal_fvq (h,den,mu,diff,q,n,pxx)
!  ==================================================
!
  integer :: h
  real    :: cal_fvq, q, n, m, d, v, fv2
  real    :: Re, Sc, X, den, mu, diff, pxx
!
  X = 0.
  if (lvent .and. h /= drop) then
    m = cal_avm(h, q, n)
    d = hydrop(h)%cm*m**(hydrop(h)%am)
    v = hydrop(h)%ctv*m**(hydrop(h)%btv)

    Re = max(min(v*d/mu, Re_max),0.)
    Sc = 0.7 !mu / diff
    X  = Sc**(1./3.) * sqrt(Re)
  endif
!
  cal_fvq = hydroc(h)%aventq + max(hydroc(h)%bventq,0.)*X + max(hydroc(h)%cventq,0.)*X*X 
!
  return
  end function
!
!  ==================================================
  function cal_fvn (h,den,mu,diff,q,n,pxx)
!  ==================================================
!
  integer :: h
  real    :: cal_fvn, q, n, m, d, v, fv2
  real    :: Re, Sc, X, den, mu, diff, pxx
!
  X = 0.
  if (lvent .and. h /= drop) then
    m = cal_avm(h, q, n)
    d = hydrop(h)%cm*m**(hydrop(h)%am)
    v = hydrop(h)%ctv*m**(hydrop(h)%btv)

    Re = max(min(v*d/mu, Re_max),0.)
    Sc = 0.7 !mu / diff
    X  = Sc**(1./3.) * sqrt(Re)
  endif
!
  cal_fvn = hydroc(h)%aventn + max(hydroc(h)%bventn,0.)*X + max(hydroc(h)%cventn,0.)*X*X
!
  return
  end function
	
end module sbpack
