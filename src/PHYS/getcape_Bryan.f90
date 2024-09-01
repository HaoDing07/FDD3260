module getcape

private

public :: get_cape

contains

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

subroutine get_cape( nk , p_in , th_in , q_in, cape, cin, z0, th0, qv0, zb, zc, src )

implicit none

integer, intent(in) :: nk
real, dimension(nk), intent(in) :: p_in,th_in,q_in
real, dimension(nk), intent(in), optional :: z0
real, intent(in), optional :: th0,qv0
integer, intent(in), optional :: src

real, intent(out) :: cape,cin
real, intent(out), optional :: zb, zc

!-----------------------------------------------------------------------
!
!  getcape - a fortran90 subroutine to calculate Convective Available
!            Potential Energy (CAPE) from a sounding.
!
!  Version 1.02                           Last modified:  10 October 2008
!
!  Author:  George H. Bryan
!           Mesoscale and Microscale Meteorology Division
!           National Center for Atmospheric Research
!           Boulder, Colorado, USA
!           gbryan@ucar.edu
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!
!  References:  Bolton (1980, MWR, p. 1046) (constants and definitions)
!               Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
!
!-----------------------------------------------------------------------
!
!  Input:     nk - number of levels in the sounding (integer)
!
!           p_in - one-dimensional array of pressure (Pa) (real)
!
!           th_in - one-dimensional array of potential temperature (K) (real)
!
!           q_in - one-dimensional array of water vapor mass mixing ratio (kg/kg) (real)
!
!  Output:  cape - Convective Available Potential Energy (J/kg) (real)
!
!            cin - Convective Inhibition (J/kg) (real)
!
!-----------------------------------------------------------------------
!  User options:

real, parameter :: pinc = 100.0   ! Pressure increment (Pa)
                              ! (smaller number yields more accurate
                              !  results,larger number makes code 
                              !  go faster)

!integer, parameter :: source = 1    ! Source parcel:
!                                ! 1 = surface
!                               ! 2 = most unstable (max theta-e)
!                                ! 3 = mixed-layer (specify ml_depth)

real, parameter :: ml_depth =  425.  ! depth (m) of mixed layer 
                                     ! for source=3

integer, parameter :: adiabat = 1   ! Formulation of moist adiabat:
                                ! 1 = pseudoadiabatic, liquid only
                                ! 2 = reversible, liquid only
                                ! 3 = pseudoadiabatic, with ice (default)
                                ! 4 = reversible, with ice

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!            No need to modify anything below here:
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

logical :: doit,ice,cloud,not_converged
integer :: k,kmax,n,nloop,i,orec,source
real, dimension(nk) :: p,t,pi,q,th,thv,z,pt,pb,pc,pn,ptv

real :: the,maxthe,parea,narea,lfc,lcl
    real :: th1,p1,t1,qv1,ql1,qi1,b1,pi1,thv1,qt,dp,dz,ps,frac
    real :: th2,p2,t2,qv2,ql2,qi2,b2,pi2,thv2,dth2
    real :: thlast,thold,dthlast,fliq,fice,tbar,qvbar,qlbar,qibar,lhv,lhs,lhf,rm,cpm
    real*8 :: avgth,avgqv

!-----------------------------------------------------------------------

    real, parameter :: g     = 9.81
    real, parameter :: p00   = 100000.0
    real, parameter :: cp    = 1004.
    real, parameter :: rd    = 288.
    real, parameter :: rv    = 438.5
    real, parameter :: xlv   = 2495800.0
    real, parameter :: xls   = 2838889.0
    real, parameter :: t0    = 273.15
    real, parameter :: cpv   = 1827.0
    real, parameter :: cpl   = 4185.0
    real, parameter :: cpi   = 2102.
    real, parameter :: lv1   = xlv !+(cpl-cpv)*t0
    real, parameter :: lv2   = cpl-cpv
    real, parameter :: ls1   = xls !+(cpi-cpv)*t0
    real, parameter :: ls2   = cpi-cpv

    real, parameter :: rp00  = 1.0/p00
    real, parameter :: eps   = 0.624052 !rd/rv
    real, parameter :: reps  = 1.6024306 !rv/rd
    real, parameter :: rddcp = rd/cp
    real, parameter :: cpdrd = cp/rd
    real, parameter :: cpdg  = cp/g

    real, parameter :: converge = 0.1

    integer, parameter :: debug_level =   0

!-----------------------------------------------------------------------

!---- Parcel source

    if (.not.present(src)) then
      source = 1
    else
      source = src
    endif

!---- convert p,t,q to mks units; get pi,q,th,thv ----!

    do k=1,nk
      p(k) = p_in(k)
      th(k) = th_in(k)
      q(k) = q_in(k) 
      pi(k) = (p(k)*rp00)**rddcp
      t(k) = th(k)*pi(k)
      thv(k) = th(k)*(1.0+reps*q(k))/(1.0+q(k))
    enddo
    if(present(zb))zb = 0.
    if(present(zc))zc = 0.

!---- get height using the hydrostatic equation ----!

    if (present(z0)) then
      z = z0
    else
      z(1) = 0.0
      do k=2,nk
        dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))
        z(k) = z(k-1) + dz
      enddo
    endif

!---- find source parcel ----!

  IF(source.eq.1)THEN
    ! use surface parcel
    kmax = 1

  ELSEIF(source.eq.2)THEN
    ! use most unstable parcel (max theta-e)

    IF(p(1).lt.50000.0)THEN
      ! first report is above 500 mb ... just use the first level reported
      kmax = 1
      maxthe = getthe(p(1),t(1),q(1))
    ELSE
      ! find max thetae below 500 mb
      maxthe = 0.0
      do k=1,nk
        if(p(k).ge.50000.0)then
          the = getthe(p(k),t(k),q(k))
          if( the.gt.maxthe )then
            maxthe = the
            kmax = k
          endif
        endif
      enddo
    ENDIF
    if(debug_level.ge.100) print *,'  kmax,maxthe = ',kmax,maxthe

  ELSEIF(source.eq.3)THEN
    ! use mixed layer

    IF( (z(2)-z(1)).gt.ml_depth )THEN
      ! the second level is above the mixed-layer depth:  just use the
      ! lowest level

      avgth = th(1)
      avgqv = q(1)
      kmax = 1

    ELSEIF( z(nk).lt.ml_depth )THEN
      ! the top-most level is within the mixed layer:  just use the
      ! upper-most level

      avgth = th(nk)
      avgqv = q(nk)
      kmax = nk

    ELSE
      ! calculate the mixed-layer properties:

      k = 2
      avgth = 0.0
      avgqv = 0.0

      do while( (z(k).le.ml_depth) .and. (k.le.nk) )
        if(debug_level.ge.100) print *,k,z(k),th(k),q(k)

        avgth = avgth + 0.5*(z(k)-z(k-1))*(th(k)+th(k-1))
        avgqv = avgqv + 0.5*(z(k)-z(k-1))*(q(k)+q(k-1))

        k = k + 1
      enddo

      th2 = th(k-1)+(th(k)-th(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))
      qv2 =  q(k-1)+( q(k)- q(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))

      avgth = avgth + 0.5*(ml_depth-z(k-1))*(th2+th(k-1))
      avgqv = avgqv + 0.5*(ml_depth-z(k-1))*(qv2+q(k-1))

      avgth = avgth/ml_depth
      avgqv = avgqv/ml_depth

      if(debug_level.ge.100) print *,k,z(k),th(k),q(k)

      kmax = 1

    ENDIF

    if(debug_level.ge.100) print *,avgth,avgqv

  ELSE

    print *
    print *,'  Unknown value for source'
    print *
    print *,'  source = ',source
    print *
    stop

  ENDIF

!---- define parcel properties at initial location ----!

    narea = 0.0
    cape = 0.0
    cin  = 0.0
    lfc  = 0.0

    k    = kmax
    pi2  = pi(kmax)
    p2   = p(kmax)
    if (present(th0).and.present(qv0)) then
      th2  = th0
      qv2  = qv0
    else if( (source.eq.1).or.(source.eq.2) )then 
      th2  = t(kmax)
      qv2  = q(kmax)
    else if( source.eq.3 )then
      th2  = avgth
      qv2  = avgqv
    endif
    t2   = th2*pi2 
    thv2 = th2*(1.0+reps*qv2)/(1.0+qv2)
    b2   = g*( thv2-thv(kmax) )/thv(kmax)

    ql2 = 0.0
    qi2 = 0.0
    qt  = qv2

    doit = .true.
    cloud = .false.
    if(adiabat.eq.1.or.adiabat.eq.2)then
      ice = .false.
    else
      ice = .true.
    endif

    the = getthe(p2,t2,qv2)
    if(debug_level.ge.100) print *,'  the = ',the

!---- begin ascent of parcel ----!

    if(debug_level.ge.100)then
      print *,'  Start loop:'
      print *,'  p2,th2,qv2 = ',p2,th2,qv2
    endif

    do while( doit .and. (k.lt.nk) )

      k = k+1
      b1 = b2

      dp = p(k-1)-p(k)

      if( dp.lt.pinc )then
        nloop = 1
      else
        nloop = 1 + int( dp/pinc )
        dp = dp/float(nloop)
      endif

      do n=1,nloop

        p1 =  p2
        t1 =  t2
        pi1 = pi2
        th1 = th2
        qv1 = qv2
        ql1 = ql2
        qi1 = qi2
        thv1 = thv2

        p2 = p2 - dp
        pi2 = (p2*rp00)**rddcp

        thold = th1
        thlast = th1
        dthlast = 0.0
        i = 0
        not_converged = .true.

        do while( not_converged )
          i = i + 1
          thold = th2
          t2 = thlast*pi2
          if(ice)then
            fliq = max(min((t2-233.15)/(273.15-233.15),1.0),0.0)
            fice = 1.0-fliq
          else
            fliq = 1.0
            fice = 0.0
          endif
          qv2 = min( qt , fliq*getqvs(p2,t2) + fice*getqvi(p2,t2) )
          qi2 = max( fice*(qt-qv2) , 0.0 )
          ql2 = max( qt-qv2-qi2 , 0.0 )

          tbar  = 0.5*(t1+t2)
          qvbar = 0.5*(qv1+qv2)
          qlbar = 0.5*(ql1+ql2)
          qibar = 0.5*(qi1+qi2)

          lhv = lv1 !- lv2*tbar
          lhs = ls1 !- ls2*tbar
          lhf = lhs - lhv

          rm=rd !+ rv*qvbar
          cpm=cp !+ cpv*qvbar + cpl*qlbar + cpi*qibar
	  
          th2=th1*exp(  lhv*(ql2-ql1)/(cpm*tbar)     &
                       +lhs*(qi2-qi1)/(cpm*tbar)     &
                       +(rm/cpm-rd/cp)*alog(p2/p1) )

          thold = th2
          dth2 = dthlast
          dthlast = th2-thlast
          if( (abs(dthlast).gt.converge).and.(abs(th2-thold).gt.converge) )then
            if( (i.gt.5).and.(sign(1.,dth2).ne.sign(1.,dthlast)) ) then
              thlast=0.5*(thlast+th2)
            else
              thlast=thlast+0.25*dthlast
            endif
            if(i.gt.45) print *,i,th2,thlast,dth2,dthlast
          else
            not_converged = .false.
          endif

          if(i.gt.50)then
            print *
            print *,'  Error:  lack of convergence'
            print *,'           Continuing        '
            print *
            exit
            !stop 1001
          endif
        enddo

        ! Latest pressure increment is complete.  Calculate some
        ! important stuff:

        if( ql2.ge.1.0e-8 ) cloud = .true.
        if(present(zc).and.zc.eq.0.and.cloud)zc = z(k)

        IF(adiabat.eq.1.or.adiabat.eq.3)THEN
          ! pseudoadiabat
          qt  = qv2
          ql2 = 0.0
          qi2 = 0.0
        ELSEIF(adiabat.le.0.or.adiabat.ge.5)THEN
          print *
          print *,'  Undefined adiabat'
          print *
          stop 10000
        ENDIF

      enddo

      thv2 = th2*(1.0+reps*qv2)/(1.0+qv2+ql2+qi2)
      b2 = g*( thv2-thv(k) )/thv(k)
      dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))
      
      the = getthe(p2,t2,qv2)

      ! Get contributions to CAPE and CIN:
      if( (b2.ge.0.0) .and. (b1.lt.0.0) )then
        ! first trip into positive area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b2/(b2-b1)
        parea =  0.5*b2*dz*frac
        narea = narea - 0.5*b1*dz*(1.0-frac)
        if(present(zb))zb = z(k)
        if(debug_level.ge.200)then
          print *,'      b1,b2 = ',b1,b2
          print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
          print *,'      frac = ',frac
          print *,'      parea = ',parea
          print *,'      narea = ',narea
        endif
        cin  = cin + narea
        narea = 0.0
      elseif( (b2.lt.0.0) .and. (b1.gt.0.0) )then
        ! first trip into neg area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b1/(b1-b2)
        parea =  0.5*b1*dz*frac
        narea = -0.5*b2*dz*(1.0-frac)
        if(debug_level.ge.200)then
          print *,'      b1,b2 = ',b1,b2
          print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
          print *,'      frac = ',frac
          print *,'      parea = ',parea
          print *,'      narea = ',narea
        endif
      elseif( b2.lt.0.0 )then
        ! still collecting negative buoyancy
        parea =  0.0
        narea = narea - 0.5*dz*(b1+b2)
      else
        ! still collecting positive buoyancy
        parea =  0.5*dz*(b1+b2)
        narea =  0.0
      endif

      cape = cape + max(0.0,parea)

      if(debug_level.ge.200)then
        write(6,102) p2,b1,b2,cape,cin,cloud
102     format(5(f13.4),2x,l1)
      endif

      if( (p(k).le.10000.0).and.(b2.lt.0.0) )then
        ! stop if b < 0 and p < 100 mb
        doit = .false.
      endif

    enddo

!---- All done ----!

    return

    contains

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    subroutine saturation_adjustment(ice,q,p,pi,t,cp,ql,qi,qv,rv,lhv)
    implicit none

    logical :: ice
    integer :: k

    real :: q,p,pi,t,ql,qi,qv
    real :: cp,rv,lhv
    real :: f,qq,qs,dq,dt

    integer, parameter :: kmax = 10
    real, parameter :: small=1.e-8

    k = 0
    dq = 1.
    qq = q
    do while ( abs(dq) > small .and. k < kmax )
     if(ice)then
        f = 1. - max(min((t-233.15)/(273.15-233.15),1.0),0.0)
      else
        f = 0.0
      endif
      qs = (1.-f)*getqvs(p,t) + f*getqvi(p,t)
      dq = (qq - qs) / ( 1. + lhv*lhv*getqvs(p,t)/(rv*cp*t*t) )
      dt = lhv*dq/cp
      qq = qq - dq
      t  = t  + dt
      k = k + 1
    enddo
    qv = qq
    ql = (1.-f)*(q - qq)
    qi = q - qq - ql

    return
    end subroutine saturation_adjustment

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real function getqvs(p,t)
    implicit none

    real :: p,t,es

    real, parameter :: eps = 0.624052

    es = 611.2*exp(17.62*(t-273.15)/(t-30.03))
    getqvs = eps*es/(p-es)

    return
    end function getqvs

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real function getqvi(p,t)
    implicit none

    real :: p,t,es

    real, parameter :: eps = 0.624052

    es = 611.2*exp(21.875*(t-273.15)/(t-7.65))
    getqvi = eps*es/(p-es)

    return
    end function getqvi

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real function getthe(p,t,q)
    implicit none

    real :: p,t,q,qs

    qs = getqvs(p,t)

    getthe=t*(100000.0/p)**(0.2868*(1.-0.24*q))*(q/qs)**(-0.4298*q*(1.-0.24*q))  &
            *exp( 2495800.*q*(1.-0.24*q)/(1004.*t) )

    !getthe=t*( (100000.0/p)**(0.2854*(1.0-0.28*q)) )   &
    !        *exp( ((3376.0/tlcl)-2.54)*q*(1.0+0.81*q) )

    return
    end function getthe

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
    end subroutine get_cape

  end module getcape
