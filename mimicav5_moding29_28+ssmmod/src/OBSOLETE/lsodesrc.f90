
#include "ctrparam.h"

! ==============================================================
!	LSODESRC:	Collection of subroutines or functions	
!			used in Chien Wang's  photochemistry
!			model from ODEPACK
!	------------------------------------------------
!  Author:
!	Chien Wang
!	MIT, Rm.E40-269, Cambridge, MA 02139
!
!  Revision:
!	Date      By			Brief Description
!	----      --			-----------------	
!	120595	  Chien Wang		created
!	070698    Chien Wang		cpp
!	072098	  Chien Wang		F90/F95
!	072098	  Chien Wang		revised
!	011700	  Chien Wang		rewote sgefa
!	070800	  Chien Wang		close xerr write for spmd runs
!	091801	  Chien Wang		removed some old style code
!	041102	  Chien Wang		do loop in f90/f95 format
!	021903	  Chien Wang		rev. for -r8
!	061604	  Chien Wang		change float() to real() remove e0
!
!	------------------------------------------------
!	Table of Contents
!
!		lsodenew:	subroutine	-- TC1
!		stodenew:	subroutine	-- TC2
!		ewset:		subroutine	-- TC3
!		solsy:		subroutine	-- TC4
!		intdy:		subroutine	-- TC5
!		cfode:		subroutine	-- TC6
!		sgefa:		subroutine	-- TC7
!		sgbfa: 		subroutine	-- TC8
!		sgesl:		subroutine	-- TC9
!		sgbsl:		subroutine	-- TC10
!		sscal:		subroutine	-- TC11
!		saxpysmp:	subroutine	-- TC12
!		xerrwv:		subroutine	-- TC13
!		r1mach:		function	-- TC14
!		vnorm:		function	-- TC15
!		isamax:		function	-- TC16
!		sdot:		function	-- TC17
!		prepj64:	subroutine	-- TC18
! ==============================================================


! ==============================================================
! --- TC1
!
      subroutine lsodenew (f, neq, y, t, tout, itol, rtol, atol,      &
                 itask,istate, iopt, rwork, lrw, iwork, liw, jac, mf)
!     ==========================================================

! ==============================================================
!   LSODENEW.F:	A simplified version of original lsode.f 
!		for cases where
! 		ISTATE = 1 &
! 		ITASK  = 1 initially
!		IOPT   = 0
!		ITOL   = 1
!        ------------------------------------------------ 
!
!  		Chien Wang
!		MIT Joint Program for Science and Policy
!			of Global Change
!
!		July 21, 1998
! ===============================================================

      external f, jac
      integer :: itol, itask, istate, iopt, lrw, liw, mf,  &
      	         iwork(liw),  neq(*) 
      real    :: t, tout, y(*), rtol(*), atol(*), rwork(lrw)

!-----------------------------------------------------------------------
! this is the march 30, 1987 version of
! lsode.. livermore solver for ordinary differential equations.
! this version is in single precision.
!
! lsode solves the initial value problem for stiff or nonstiff
! systems of first order ode-s,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
! lsode is a package based on the gear and gearb packages, and on the
! october 23, 1978 version of the tentative odepack user interface
! standard, with minor modifications.
!-----------------------------------------------------------------------
! reference..
!     alan c. hindmarsh,  odepack, a systematized collection of ode
!     solvers, in scientific computing, r. s. stepleman et al. (eds.),
!     north-holland, amsterdam, 1983, pp. 55-64.
!-----------------------------------------------------------------------
! author and contact.. alan c. hindmarsh,
!                      computing and mathematics research div., l-316
!                      lawrence livermore national laboratory
!                      livermore, ca 94550.
!-----------------------------------------------------------------------

      external prepj, solsy
      integer :: iii, illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,        &
     		mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns,	 &
     		icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,  &
     		maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu,  &
     		i, i1, i2, iflag, imxer, kgo, lf0,			 &
     		leniw, lenrw, lenwm, ml, mu, mxhnl0, mxstp0,		 &
     		mord(2) 						  
     real    :: rowns, ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,	 &
     		atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,	 &
     		tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0,	 &
     		r1mach, vnorm, r_tmp

      logical :: ihit
!-----------------------------------------------------------------------
! the following internal common block contains
! (a) variables which are local to any subroutine but whose values must
!     be preserved between calls to the routine (own variables), and
! (b) variables which are communicated between subroutines.
! the structure of the block is as follows..  all real variables are
! listed first, followed by all integers.  within each type, the
! variables are grouped with those local to subroutine lsode first,
! then those local to subroutine stode, and finally those used
! for communication.  the block is declared in subroutines
! lsode, intdy, stode, prepj, and solsy.  groups of variables are
! replaced by dummy arrays in the common declarations in routines
! where those variables are not used.
!-----------------------------------------------------------------------
	common /ls0001/ rowns(209),                               &
     	ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,		  &
     	illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,	  &
     	mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),	  &
     	icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
     	maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
!
!      data  mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
      data  mord(1),mord(2)/12,5/, mxstp0/260/, mxhnl0/10/
      data illin/0/, ntrep/0/
!-----------------------------------------------------------------------
! block a.
! this code block is executed on every call.
! it tests istate and itask for legality and branches appropriately.
! if istate .gt. 1 but the flag init shows that initialization has
! not yet been done, an error return occurs.
! if istate = 1 and tout = t, jump to block g and return immediately.
!-----------------------------------------------------------------------

      init = 0
      ntrep = 0
!------------------------------------------------------------------------
! block b.
! the next code block is executed for the initial call (istate = 1),
! or for a continuation call with parameter changes (istate = 3).
! it contains checking of all inputs and various initializations.
!
! first check legality of the non-optional inputs neq, itol, iopt,
! mf, ml, and mu.
!-----------------------------------------------------------------------

      n = neq(1)
      meth = mf/10
      miter = mf - 10*meth
      if (miter .gt. 3)then
        ml = iwork(1)
        mu = iwork(2)
      endif

! next process and check the optional inputs. --------------------------
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      h0 = 0.0
      hmxi = 0.0
      hmin = 0.0

!-----------------------------------------------------------------------
! set work array pointers and check lengths lrw and liw.
! pointers to segments of rwork and iwork are named by prefixing l to
! the name of the segment.  e.g., the segment yh starts at rwork(lyh).
! segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.
!-----------------------------------------------------------------------

	lyh = 21
	if (istate .eq. 1) nyh = n
	lwm = lyh + (maxord + 1)*nyh
	select case (miter)
      	  case (0)
		lenwm = 0
	  case (1, 2)
		lenwm = n*n + 2
	  case (3)
		lenwm = n + 2
	  case (4)
		lenwm = (2*ml + mu + 1)*n + 2
	end select
	lewt = lwm + lenwm
	lsavf = lewt + n
	lacor = lsavf + n
	lenrw = lacor + n - 1
	iwork(17) = lenrw
	liwm = 1
	leniw = 20 + n
	if (miter .eq. 0 .or. miter .eq. 3) leniw = 20
	iwork(18) = leniw

! check rtol and atol for legality. ------------------------------------
      rtoli = rtol(1)
      atoli = atol(1)

!      if(itol.eq.2)then
!      do 70 i = 1,n
!        atoli = atol(i)
!70    continue
!      endif

!-----------------------------------------------------------------------
! block c.
! the next block is for the initial call only (istate = 1).
! it contains all remaining initializations, the initial call to f,
! and the calculation of the initial step size.
! the error weights in ewt are inverted after being loaded.
!-----------------------------------------------------------------------

      uround = r1mach(4)
      !uround = 1.e-20
      tn = t
      jstart = 0
      if (miter .gt. 0) rwork(lwm) = sqrt(uround)
      nhnil = 0
      nst = 0
      nje = 0
      nslast = 0
      hu = 0.0
      nqu = 0
      ccmax = 0.3
      maxcor = 3
      msbp = 20
      mxncf = 10

! initial call to f.  (lf0 points to yh(*,2).) -------------------------
      lf0 = lyh + nyh
      call f (neq, t, y, rwork(lf0))
      nfe = 1

! load the initial value vector in yh. ---------------------------------
      iii = lyh-1
      do i = 1,n
	rwork(i+iii) = y(i)
      end do

! load and invert the ewt array.  (h is temporarily set to 1.0.) -------
      nq = 1
      h = 1.0
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))

      iii = lewt-1
      do i = lewt,iii+n
	rwork(i) = 1.0/rwork(i)
      end do
!-----------------------------------------------------------------------
! the coding below computes the step size, h0, to be attempted on the
! first step, unless the user has supplied a value for this.
! first check that tout - t differs significantly from zero.
! a scalar tolerance quantity tol is computed, as max(rtol(i))
! if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
! so as to be between 100*uround and 1.0e-3.
! then the computed value h0 is given by..
!                                      neq
!   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
!                                       1
! where   w0     = max ( abs(t), abs(tout) ),
!         f(i)   = i-th component of initial value of f,
!         ywt(i) = ewt(i)/tol  (a weight for y(i)).
! the sign of h0 is inferred from the initial values of tout and t.
!-----------------------------------------------------------------------

      if (h0 .eq. 0.0)then

        tdist = abs(tout - t)
        w0 = max(abs(t),abs(tout))
        tol = rtol(1)

!        if (itol .gt. 2)then
!      do 130 i = 1,n
! 130    tol = max(tol,rtol(i))
!        endif

	if (tol .le. 0.0)then
          atoli = atol(1)
	do i = 1,n
!        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
	  ayi = abs(y(i))
	  if (ayi .ne. 0.0) tol = max(tol,atoli/ayi)
	end do
        endif

	tol = max(tol,100.0*uround)
	tol = min(tol,0.001)
	sum = vnorm (n, rwork(lf0), rwork(lewt))
	sum = 1.0/(tol*w0*w0) + tol*sum**2
	h0 = 1.0/sqrt(sum)
	h0 = min(h0,tdist)
	h0 = sign(h0,tout-t)

      endif

! adjust h0 if necessary to meet hmax bound. ---------------------------
      rh = abs(h0)*hmxi
      if (rh .gt. 1.0) h0 = h0/rh
! load h with h0 and scale yh(*,2) by h0. ------------------------------
      h = h0
      do i = 1,n
	rwork(i+lf0-1) = h0*rwork(i+lf0-1)
      end do

!-----------------------------------------------------------------------
! block e.
! the next block is normally executed for all calls and contains
! the call to the one-step core integrator stode.
!
! this is a looping point for the integration steps.
!
! first check for too many steps being taken, update ewt (if not at
! start of problem), check for too much accuracy being requested, and
! check for h below the roundoff level in t.
!-----------------------------------------------------------------------
      
      do 270 while ((tn - tout)*h .lt. 0.0)

      tolsf = uround*vnorm (n, rwork(lyh), rwork(lewt))

      if (tolsf .gt. 1.0)then
       tolsf = tolsf*2.0
       go to 580
      endif

      call stodenew (neq, y, rwork(lyh), nyh,                  &
     		    rwork(lyh), rwork(lewt),		       &
     	rwork(lsavf), rwork(lacor), rwork(lwm), iwork(liwm),   &
     	f, jac, prepj, solsy)

      init = 1

      if ((nst-nslast) .ge. mxstep) go to 580

      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))

      do 260 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0) go to 580
	rwork(i+lewt-1) = 1.0/rwork(i+lewt-1)
260   continue

270   continue

      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout

      istate = 2
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return


! --- set y vector, t, illin, and optional outputs.
580	do i = 1,n
	  y(i) = rwork(i+lyh-1)
	end do
	t = tn
	illin = 0
	rwork(11) = hu
	rwork(12) = h
	rwork(13) = tn
	iwork(11) = nst
	iwork(12) = nfe
	iwork(13) = nje
	iwork(14) = nqu
	iwork(15) = nq
	
	return
	 end


! ==============================================================
! -- TC2
!
      subroutine stodenew (neq, y, yh, nyh, yh1, ewt, savf, acor,   &
        wm, iwm, f, jac, pjac, slvs)
!     ==========================================================

!---------------------------------------------------------------
!  STODENEW.F:	A simplified version of STODE.F	
!			for JSTART >= 0	
!	--------------------------------------------
!
!  		Chien Wang
!
!  Last revised:	July 21, 1998
!
!---------------------------------------------------------------

!lll. optimize
      external f, jac, pjac, slvs
      integer :: neq(*), nyh, iwm(*),                                  &
     		iownd, ialth, ipup, lmax, meo, nqnyh, nslp,	       &
     		icf, ierpj, iersl, jcur, jstart, kflag, l, 	       &
     		meth, miter, maxord, maxcor, msbp,		       &
     		mxncf, n, nq, nst, nfe, nje, nqu,		       &
     		i, i1, iredo, iret, j, jb, m, ncf, newq 	        
     real    :: y(*), yh(nyh,*), yh1(*), ewt(*), savf(*), acor(*),     &
     		wm(*), conit, crate, el, elco, hold, rmax, tesco,      &
     		ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,	       &
     		dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,      &
     		r, rh, rhdn, rhsm, rhup, told, vnorm		        
     common /ls0001/ conit, crate, el(13), elco(13,12), 	       &
     	hold, rmax, tesco(3,12),				       &
     	ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14),      &
     	ialth, ipup, lmax, meo, nqnyh, nslp,			       &
     	icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,        &
     	maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
! -----------------------------------------------------------------------
! stode performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! note.. stode is independent of the value of the iteration method
! indicator miter, when this is .ne. 0, and hence is independent
! of the type of chord method used, or the jacobian structure.
! communication with stode is done with the following variables..
!
! neq    = integer array containing problem size in neq(1), and
!          passed as the neq argument in all calls to f and jac.
! y      = an array of length .ge. n used as the y argument in
!          all calls to f and jac.
! yh     = an nyh by lmax array containing the dependent variables
!          and their approximate scaled derivatives, where
!          lmax = maxord + 1.  yh(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by h**j/factorial(j)
!          (j = 0,1,...,nq).  on entry for the first step, the first
!          two columns of yh must be set from the initial values.
! nyh    = a constant integer .ge. n, the first dimension of yh.
! yh1    = a one-dimensional array occupying the same space as yh.
! ewt    = an array of length n containing multiplicative weights
!          for local error measurements.  local errors in y(i) are
!          compared to 1.0/ewt(i) in various error tests.
! savf   = an array of working storage, of length n.
!          also used for input of yh(*,maxord+2) when jstart = -1
!          and maxord .lt. the current order nq.
! acor   = a work array of length n, used for the accumulated
!          corrections.  on a successful return, acor(i) contains
!          the estimated one-step local error in y(i).
! wm,iwm = real and integer work arrays associated with matrix
!          operations in chord iteration (miter .ne. 0).
! pjac   = name of routine to evaluate and preprocess jacobian matrix
!          and p = i - h*el0*jac, if a chord method is being used.
! slvs   = name of routine to solve linear system in chord iteration.
! ccmax  = maximum relative change in h*el0 before pjac is called.
! h      = the step size to be attempted on the next step.
!          h is altered by the error control algorithm during the
!          problem.  h can be either positive or negative, but its
!          sign must remain constant throughout the problem.
! hmin   = the minimum absolute value of the step size h to be used.
! hmxi   = inverse of the maximum absolute value of h to be used.
!          hmxi = 0.0 is allowed and corresponds to an infinite hmax.
!          hmin and hmxi may be changed at any time, but will not
!          take effect until the next change of h is considered.
! tn     = the independent variable. tn is updated on each step taken.
! jstart = an integer used for input only, with the following
!          values and meanings..
!               0  perform the first step.
!           .gt.0  take a new step continuing from the last.
!              -1  take the next step with a new value of h, maxord,
!                    n, meth, miter, and/or matrix parameters.
!              -2  take the next step with a new value of h,
!                    but with other inputs unchanged.
!          on return, jstart is set to 1 to facilitate continuation.
! kflag  = a completion code with the following meanings..
!               0  the step was succesful.
!              -1  the requested error could not be achieved.
!              -2  corrector convergence could not be achieved.
!              -3  fatal error in pjac or slvs.
!          a return with kflag = -1 or -2 means either
!          abs(h) = hmin or 10 consecutive failures occurred.
!          on a return with kflag negative, the values of tn and
!          the yh array are as of the beginning of the last
!          step, and h is the last step size attempted.
! maxord = the maximum order of integration method to be allowed.
! maxcor = the maximum number of corrector iterations allowed.
! msbp   = maximum number of steps between pjac calls (miter .gt. 0).
! mxncf  = maximum number of convergence failures allowed.
! meth/miter = the method flags.  see description in driver.
! n      = the number of first-order differential equations.
! -----------------------------------------------------------------------
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0
! -----------------------------------------------------------------------
! on the first call, the order is set to 1, and other variables are
! initialized.  rmax is the maximum ratio by which h can be increased
! in a single step.  it is initially 1.e4 to compensate for the small
! initial h, but then is normally equal to 10.  if a failure
! occurs (in corrector convergence or error test), rmax is set at 2
! for the next increase.
!
! cfode is called to get all the integration coefficients for the
! current meth.  then the el vector and related constants are reset
! whenever the order nq is changed, or at the start of the problem.
! -----------------------------------------------------------------------
      if(jstart.eq.0)then
        lmax = maxord + 1
        nq = 1
        l = 2
        ialth = 2
        rmax = 10000.0
        rc = 0.0
        el0 = 1.0
        crate = 0.7
        hold = h
        meo = meth
        nslp = 0
        ipup = miter
        iret = 3

      	call cfode (meth, elco, tesco)
      endif
! -----------------------------------------------------------------------

      if (jstart .gt. 0) go to 200

 150  do i = 1,l
	el(i) = elco(i,nq)
      end do
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5/real(nq+2)

      if(iret.eq.3) go to 200

! -----------------------------------------------------------------------
! if h is being changed, the h ratio rh is checked against
! rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to
! l = nq + 1 to prevent a change of h for that many steps, unless
! forced by a convergence or error test failure.
! -----------------------------------------------------------------------

 170  rh = max(rh,hmin/abs(h))
      rh = min(rh,rmax)
      rh = rh/max(1.0,abs(h)*hmxi*rh)
      r = 1.0
      do j = 2,l
        r = r*rh
        do i = 1,n
	  yh(i,j) = yh(i,j)*r
	end do
      end do      
      h = h*rh
      rc = rc*rh
      ialth = l
      if (iredo .eq. 0) go to 690
! -----------------------------------------------------------------------
! this section computes the predicted values by effectively
! multiplying the yh array by the pascal triangle matrix.
! rc is the ratio of new to old values of the coefficient  h*el(1).
! when rc differs from 1 by more than ccmax, ipup is set to miter
! to force pjac to be called, if a jacobian is involved.
! in any case, pjac is called at least every msbp steps.
! -----------------------------------------------------------------------
 200  if (abs(rc-1.0) .gt. ccmax) ipup = miter
      if (nst .ge. nslp+msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      do jb = 1,nq
        i1 = i1 - nyh
!dir$ ivdep
        do i = i1,nqnyh
	  yh1(i) = yh1(i) + yh1(i+nyh)
 	end do
      end do
! -----------------------------------------------------------------------
! up to maxcor corrector iterations are taken.  a convergence test is
! made on the r.m.s. norm of each correction, weighted by the error
! weight vector ewt.  the sum of the corrections is accumulated in the
! vector acor(i).  the yh array is not altered in the corrector loop.
! -----------------------------------------------------------------------
220	m = 0
	y(1:n) = yh(1:n,1)
	call f (neq, tn, y, savf)
	nfe = nfe + 1
! -----------------------------------------------------------------------
! if indicated, the matrix p = i - h*el(1)*j is reevaluated and
! preprocessed before starting the corrector iteration.  ipup is set
! to 0 as an indicator that this has been done.
! -----------------------------------------------------------------------
	if (ipup .gt. 0)then
	  call pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)
	  ipup = 0
	  rc = 1.0
	  nslp = nst
	  crate = 0.7
	  if (ierpj .ne. 0) go to 430
	endif
	
	acor(1:n) = 0.0
	
270	if (miter .ne. 0) go to 350
! -----------------------------------------------------------------------
! in the case of functional iteration, update y directly from
! the result of the last function evaluation.
! -----------------------------------------------------------------------
	do i = 1,n
	  savf(i) = h*savf(i) - yh(i,2)
	  y(i)    = savf(i) - acor(i)
	end do
	del = vnorm (n, y, ewt)
	do i = 1,n
	  y(i)    = yh(i,1) + el(1)*savf(i)
	  acor(i) = savf(i)
	end do
	go to 400
! -----------------------------------------------------------------------
! in the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! p as coefficient matrix.
! -----------------------------------------------------------------------
350	y(1:n) = h*savf(1:n) - (yh(1:n,2) + acor(1:n))
	call slvs (wm, iwm, y, savf)
	if (iersl .lt. 0) go to 430
	if (iersl .gt. 0) go to 410
	del = vnorm (n, y, ewt)
	do i = 1,n
	  acor(i) = acor(i) + y(i)
	  y(i) = yh(i,1) + el(1)*acor(i)
	end do
! -----------------------------------------------------------------------
! test for convergence.  if m.gt.0, an estimate of the convergence
! rate constant is stored in crate, and this is used in the test.
! -----------------------------------------------------------------------
 400  if (m .ne. 0) crate = max(0.2*crate,del/delp)
      dcon = del*min(1.0,1.5*crate)/(tesco(2,nq)*conit)
      if (dcon .le. 1.0) go to 450
      m = m + 1
      if (m .eq. maxcor) go to 410
      if (m .ge. 2 .and. del .gt. 2.0*delp) go to 410
      delp = del
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      go to 270
! -----------------------------------------------------------------------
! the corrector iteration failed to converge.
! if miter .ne. 0 and the jacobian is out of date, pjac is called for
! the next try.  otherwise the yh array is retracted to its values
! before prediction, and h is reduced, if possible.  if h cannot be
! reduced or mxncf failures have occurred, exit with kflag = -2.
! -----------------------------------------------------------------------
 410  if (miter .eq. 0 .or. jcur .eq. 1) go to 430
      icf = 1
      ipup = miter
      go to 220
 430  icf = 2
      ncf = ncf + 1
      rmax = 2.0
      tn = told
      i1 = nqnyh + 1
	do jb = 1,nq
	  i1 = i1 - nyh
!dir$ ivdep
	  do i = i1,nqnyh
		yh1(i) = yh1(i) - yh1(i+nyh)
	  end do
	end do
      if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
      if (abs(h) .le. hmin*1.00001) go to 670
      if (ncf .eq. mxncf) go to 670
      rh = 0.25
      ipup = miter
      iredo = 1
      go to 170
! -----------------------------------------------------------------------
! the corrector has converged.  jcur is set to 0
! to signal that the jacobian involved may need updating later.
! the local error test is made and control passes to statement 500
! if it fails.
! -----------------------------------------------------------------------
 450  jcur = 0
      if (m .eq. 0) dsm = del/tesco(2,nq)
      if (m .gt. 0) dsm = vnorm (n, acor, ewt)/tesco(2,nq)
      if (dsm .gt. 1.0) go to 500
! -----------------------------------------------------------------------
! after a successful step, update the yh array.
! consider changing h if ialth = 1.  otherwise decrease ialth by 1.
! if ialth is then 1 and nq .lt. maxord, then acor is saved for
! use in a possible order increase on the next step.
! if a change in h is considered, an increase or decrease in order
! by one is considered also.  a change in h is made only if it is by a
! factor of at least 1.1.  if not, ialth is set to 3 to prevent
! testing for that many steps.
! -----------------------------------------------------------------------
	kflag = 0
	iredo = 0
	nst = nst + 1
	hu = h
	nqu = nq
	do j = 1,l
        do i = 1,n
	  yh(i,j) = yh(i,j) + el(j)*acor(i)
	end do
	end do
	ialth = ialth - 1
	if (ialth .eq. 0) go to 520
	if (ialth .gt. 1) go to 700
	if (l .eq. lmax) go to 700
	yh(1:n,lmax) = acor(1:n)
	go to 700
! -----------------------------------------------------------------------
! the error test failed.  kflag keeps track of multiple failures.
! restore tn and the yh array to their previous values, and prepare
! to try the step again.  compute the optimum step size for this or
! one lower order.  after 2 or more failures, h is forced to decrease
! by a factor of 0.2 or less.
! -----------------------------------------------------------------------
 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      do 515 jb = 1,nq
        i1 = i1 - nyh
!dir$ ivdep
        do 510 i = i1,nqnyh
	  yh1(i) = yh1(i) - yh1(i+nyh)
 510	continue
 515    continue
      rmax = 2.0
      if (abs(h) .le. hmin*1.00001) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      rhup = 0.0
      go to 540
! -----------------------------------------------------------------------
! regardless of the success or failure of the step, factors
! rhdn, rhsm, and rhup are computed, by which h could be multiplied
! at order nq - 1, order nq, or order nq + 1, respectively.
! in the case of failure, rhup = 0.0 to avoid an order increase.
! the largest of these is determined and the new order chosen
! accordingly.  if the order is to be increased, we compute one
! additional scaled derivative.
! -----------------------------------------------------------------------
520	rhup = 0.0
	if (l .eq. lmax) go to 540
	savf(1:n) = acor(1:n) - yh(1:n,lmax)
	dup = vnorm (n, savf, ewt)/tesco(3,nq)
	exup = 1.0/real(l+1)
	rhup = 1.0/(1.4*dup**exup + 0.0000014)
540	exsm = 1.0/real(l)
	rhsm = 1.0/(1.2*dsm**exsm + 0.0000012)
	rhdn = 0.0
	if (nq .eq. 1) go to 560
	ddn = vnorm (n, yh(1,l), ewt)/tesco(1,nq)
	exdn = 1.0/real(nq)
	rhdn = 1.0/(1.3*ddn**exdn + 0.0000013)
560	if (rhsm .ge. rhup) go to 570
	if (rhup .gt. rhdn) go to 590
	go to 580
570	if (rhsm .lt. rhdn) go to 580
	newq = nq
	rh = rhsm
	go to 620
580	newq = nq - 1
	rh = rhdn
	if (kflag .lt. 0 .and. rh .gt. 1.0) rh = 1.0
	go to 620
590	newq = l
	rh = rhup
	if (rh .lt. 1.1) go to 610
	r = el(l)/real(l)
	yh(1:n,newq+1) = acor(1:n)*r
	go to 630
610	ialth = 3
	go to 700
620	if ((kflag .eq. 0) .and. (rh .lt. 1.1)) go to 610
	if (kflag .le. -2) rh = min(rh,0.2)
! -----------------------------------------------------------------------
! if there is a change of order, reset nq, l, and the coefficients.
! in any case h is reset according to rh and the yh array is rescaled.
! then exit from 690 if the step was ok, or redo the step otherwise.
! -----------------------------------------------------------------------
      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150
! -----------------------------------------------------------------------
! control reaches this section if 3 or more failures have occured.
! if 10 failures have occurred, exit with kflag = -1.
! it is assumed that the derivatives that have accumulated in the
! yh array have errors of the wrong order.  hence the first
! derivative is recomputed, and the order is set to 1.  then
! h is reduced by a factor of 10, and the step is retried,
! until it succeeds or h reaches hmin.
! -----------------------------------------------------------------------
 640  if (kflag .eq. -10) go to 660
      rh = 0.1
      rh = max(hmin/abs(h),rh)
      h = h*rh
      do 645 i = 1,n
	y(i) = yh(i,1)
 645  continue
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      do 650 i = 1,n
        yh(i,2) = h*savf(i)
 650  continue
      ipup = miter
      ialth = 5
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150
! -----------------------------------------------------------------------
! all returns are made through this section.  h is saved in hold
! to allow the caller to change h on the next step.
! -----------------------------------------------------------------------
660	kflag = -1
	go to 720
670	kflag = -2
	go to 720
680	kflag = -3
	go to 720
690	rmax = 10.0
700	r = 1.0/tesco(2,nqu)
	acor(1:n) = acor(1:n)*r
720	hold = h
	jstart = 1
	
	return
	 end

! =====================================================
! -- TC3
!
      subroutine ewset (n, itol, rtol, atol, ycur, ewt)
!     =================================================

!lll. optimize
! --------------------------------------------------------------------
! this subroutine sets the error weight vector ewt according to
!     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
! with the subscript on rtol and/or atol possibly replaced by 1 above,
! depending on the value of itol.
! --------------------------------------------------------------------
!
! --- Chien Wang 
! --- rewritten 031795, 072198
!
	integer :: n, itol, i
	real    :: rtol(*), atol(*), ycur(n), ewt(n)

	ewt = rtol(1)*abs(ycur) + atol(1)

	return
	 end

! ======================================
! -- TC4
!
      subroutine solsy (wm, iwm, x, tem)
!     =================================

!
! --- Chien Wang
! --- MIT
! --- Rewrote 122695, 072198
!

!lll. optimize
	integer :: iwm(*), iownd, iowns, icf, ierpj, iersl,           &
     		  jcur, jstart, kflag, l, meth, miter,		      &
     		  maxord, maxcor, msbp, mxncf, n, nq, nst, 	      &
     		  nfe, nje, nqu, i, meband, ml, mu
       real    :: wm(*), x(*), tem(*), rowns, ccmax, el0, h,	      &
     		  hmin, hmxi, hu, rc, tn, uround,		      &
     		  di, hl0, phl0, r

       common /ls0001/ rowns(209),				      &
     	ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,		      &
     	iownd(14), iowns(6),					      &
     	icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,       &
     	maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
! ---------------------------------------------------------------------
! this routine manages the solution of the linear system arising from
! a chord iteration.  it is called if miter .ne. 0.
! if miter is 1 or 2, it calls sgesl to accomplish this.
! if miter = 3 it updates the coefficient h*el0 in the diagonal
! matrix, and then computes the solution.
! if miter is 4 or 5, it calls sgbsl.
! communication with solsy uses the following variables..
! wm    = real work space containing the inverse diagonal matrix if
!         miter = 3 and the lu decomposition of the matrix otherwise.
!         storage of matrix elements starts at wm(3).
!         wm also contains the following matrix-related data..
!         wm(1) = sqrt(uround) (not used here),
!         wm(2) = hl0, the previous value of h*el0, used if miter = 3.
! iwm   = integer work space containing pivot information, starting at
!         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band
!         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
! x     = the right-hand side vector on input, and the solution vector
!         on output, of length n.
! tem   = vector of work space of length n, not used in this version.
! iersl = output flag (in common).  iersl = 0 if no trouble occurred.
!         iersl = 1 if a singular matrix arose with miter = 3.
! this routine also uses the common variables el0, h, miter, and n.
! ---------------------------------------------------------------------
	
	iersl = 0
	
	select case(miter)
	  case (1, 2)
		call sgesl (wm(3), n, n, iwm(21), x, 0)
	  case (3)
	  	call solsy2(wm, iwm, x, tem)
	  case (4, 5)
	  	ml = iwm(1)
		mu = iwm(2)
		meband = 2*ml + mu + 1
		call sgbsl(wm(3),meband,n,ml,mu,iwm(21),x,0)
	end select 
	
	return
	 end

!	==================================
	subroutine solsy2(wm, iwm, x, tem)
!	==================================

! --- Chien Wang
! --- MIT
! --- 122695, 072198

	integer :: iwm(*), iownd, iowns, icf, ierpj, iersl,    &
     		  jcur, jstart, kflag, l, meth, miter,	       &
     		  maxord, maxcor, msbp, mxncf, n, nq, nst,     &
     		  nfe, nje, nqu, i, meband, ml, mu
       real    :: wm(*), x(*), tem(*), rowns, ccmax, el0, h,   &
     		  hmin, hmxi, hu, rc, tn, uround,	       &
     		  di, hl0, phl0, r

       common /ls0001/ rowns(209),			       &
     	ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,	       &
     	iownd(14), iowns(6),				       &
     	icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,&
     	maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu

	iersl = 0
	phl0 = wm(2)
	hl0 = h*el0
	wm(2) = hl0

	if (hl0 .ne. phl0)then
	  r = hl0/phl0
	  do i = 1,n
		di = 1.0 - r*(1.0 - 1.0/wm(i+2))
		if (abs(di) .ne. 0.0)then
		  wm(i+2) = 1.0/di
		else
	          iersl = 1
	          goto 360
	   	endif
	  end do
	else
	  do i = 1,n
       	      x(i) = wm(i+2)*x(i)
	  end do
	endif

360     return
      	 end

! ================================================
! -- TC5
!
      subroutine intdy (t, k, yh, nyh, dky, iflag)
!     ===========================================

!lll. optimize
      integer k, nyh, iflag
      integer iownd, iowns,                                       &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,  &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      real t, yh, dky
      real rowns,ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real c, r, s, tp
      dimension yh(nyh,*), dky(*)
      common /ls0001/ rowns(209),                                &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,          &
         iownd(14), iowns(6),                                    &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
! ----------------------------------------------------------------------
! intdy computes interpolated values of the k-th derivative of the
! dependent variable vector y, and stores it in dky.  this routine
! is called within the package with k = 0 and t = tout, but may
! also be called by the user for any k up to the current order.
! (see detailed instructions in the usage documentation.)
! ----------------------------------------------------------------------
! the computed values in dky are gotten by interpolation using the
! nordsieck history array yh.  this array corresponds uniquely to a
! vector-valued polynomial of degree nqcur or less, and dky is set
! to the k-th derivative of this polynomial at t.
! the formula for dky is..
!              q
!  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
!             j=k
! where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
! the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are
! communicated by common.  the above sum is done in reverse order.
! iflag is returned negative if either k or t is out of bounds.
! ----------------------------------------------------------------------
      iflag = 0
      if (k .lt. 0 .or. k .gt. nq) go to 80
      tp = tn - hu -  100.0*uround*(tn + hu)
      if ((t-tp)*(t-tn) .gt. 0.0) go to 90
!
      s = (t - tn)/h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = l - k
      do 10 jj = jj1,nq
        ic = ic*jj
 10   continue
 15   c = real(ic)
      do 20 i = 1,n
        dky(i) = c*yh(i,l)
 20   continue
      if (k .eq. nq) go to 55
      jb2 = nq - k
      do 50 jb = 1,jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        if (k .eq. 0) go to 35
        jj1 = jp1 - k
        do 30 jj = jj1,j
          ic = ic*jj
 30     continue
 35     c = real(ic)
        do 40 i = 1,n
          dky(i) = c*yh(i,jp1) + s*dky(i)
 40	continue
 50     continue
      if (k .eq. 0) return
 55   r = h**(-k)
      do 60 i = 1,n
        dky(i) = r*dky(i)
 60   continue
      return
!
 80   call xerrwv(30hintdy--  k (=i1) illegal      ,  &
        30, 51, 0, 1, k, 0, 0, 0.0, 0.0)
      iflag = -1
      return
 90   call xerrwv(30hintdy--  t (=r1) illegal      , &
        30, 52, 0, 0, 0, 0, 1, t, 0.0)
      call xerrwv(    &
       60h      t not in interval tcur - hu (= r1) to tcur (=r2)      ,  &
        60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      
	return
	 end

! ========================================
! -- TC6
!
      subroutine cfode (meth, elco, tesco)
!     ===================================

! --- Chien Wang
! --- MIT
! --- rewritten 072198

!lll. optimize
	integer :: meth, i, ib, nq, nqm1, nqp1
	real    :: elco(13,12), tesco(3,12),         &
     		  agamq, fnq, fnqm1, pc, pint, ragq, &
     		  rqfac, rq1fac, tsign, xpin

! ---------------------------------------------------------------------
! cfode is called by the integrator routine to set coefficients
! needed there.  the coefficients for the current method, as
! given by the value of meth, are set for all orders and saved.
! the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.
! (a smaller value of the maximum order is also allowed.)
! cfode is called once at the beginning of the problem,
! and is not called again unless and until meth is changed.
!
! the elco array contains the basic method coefficients.
! the coefficients el(i), 1 .le. i .le. nq+1, for the method of
! order nq are stored in elco(i,nq).  they are given by a genetrating
! polynomial, i.e.,
!     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
! for the implicit adams methods, l(x) is given by
!     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
! for the bdf methods, l(x) is given by
!     l(x) = (x+1)*(x+2)* ... *(x+nq)/k,
! where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
!
! the tesco array contains test constants used for the
! local error test and the selection of step size and/or order.
! at order nq, tesco(k,nq) is used for the selection of step
! size at order nq - 1 if k = 1, at order nq if k = 2, and at order
! nq + 1 if k = 3.
! ---------------------------------------------------------------------

	select case (meth)
	  case (1)
	  	call cfode1(meth, elco, tesco)
	  case (2)
	  	call cfode2(meth, elco, tesco)
	end select

	return
	 end

!	==================================
	subroutine cfode1(meth,elco,tesco)
!	==================================


! --- Chien Wang
! --- MIT
! --- 122695, 072198

	integer :: meth, i, ib, nq, nqm1, nqp1
	real    :: elco(13,12), tesco(3,12),               &
     		  agamq, fnq, fnqm1, pc(12), pint, ragq,   &
     		  rqfac, rq1fac, tsign, xpin

	elco (1:2,1) = 1.0
	tesco(1,1)   = 0.0
	tesco(2,1)   = 2.0
	tesco(1,2)   = 1.0
	tesco(3,12)  = 0.0
	pc(1)        = 1.0
	rqfac        = 1.0
	
	do 140 nq = 2,12
! -------------------------------------------------------------
! the pc array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq-1).
! initially, p(x) = 1.
! -------------------------------------------------------------
	rq1fac = rqfac
        rqfac = rqfac/real(nq)
        nqm1 = nq - 1
        fnqm1 = real(nqm1)
        nqp1 = nq + 1
! --- form coefficients of p(x)*(x+nq-1).
        pc(nq) = 0.0
        do ib = 1,nqm1
          i = nqp1 - ib
	  pc(i) = pc(i-1) + fnqm1*pc(i)
 	end do
        pc(1) = fnqm1*pc(1)
! --- compute integral, -1 to 0, of p(x) and x*p(x).
        pint = pc(1)
        xpin = pc(1)/2.0
        tsign = 1.0
        do i = 2,nq
          tsign = -tsign
          pint = pint + tsign*pc(i)/real(i)
	  xpin = xpin + tsign*pc(i)/real(i+1)
	end do
! --- store coefficients in elco and tesco.
        elco(1,nq) = pint*rq1fac
        elco(2,nq) = 1.0
        do i = 2,nq
	  elco(i+1,nq) = rq1fac*pc(i)/real(i)
 	end do
        agamq = rqfac*xpin
        ragq = 1.0/agamq
        tesco(2,nq) = ragq
        if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/real(nqp1)
        tesco(3,nqm1) = ragq
 140    continue
 
	return
	 end

!	====================================
	subroutine cfode2(meth, elco, tesco)
!	====================================

! --- Chien Wang
! --- MIT
! --- 122695

	integer :: meth, i, ib, nq, nqm1, nqp1 
	real    :: elco(13,12), tesco(3,12),               &
     		  agamq, fnq, fnqm1, pc(12), pint, ragq,   &
     		  rqfac, rq1fac, tsign, xpin
 
	pc(1) = 1.0
	rq1fac = 1.0
	
	do 230 nq = 1,5
! ------------------------------------------------------------
! the pc array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq).
! initially, p(x) = 1.
! ------------------------------------------------------------
        fnq = real(nq)
        nqp1 = nq + 1
! --- form coefficients of p(x)*(x+nq).
        pc(nqp1) = 0.0
        do ib = 1,nq
          i = nq + 2 - ib
	  pc(i) = pc(i-1) + fnq*pc(i)
 	end do
        pc(1) = fnq*pc(1)
! --- store coefficients in elco and tesco.
	elco(1:nqp1,nq) = pc(1:nqp1)/pc(2)
        elco(2,nq) = 1.0
        tesco(1,nq) = rq1fac
        tesco(2,nq) = real(nqp1)/elco(1,nq)
        tesco(3,nq) = real(nq+2)/elco(1,nq)
        rq1fac = rq1fac/fnq
 230    continue
	
	return
	 end

! =======================================
! -- TC7
!
      subroutine sgefa(a,lda,n,ipvt,info)
!     ==================================

! --- Chien Wang
! --- MIT
! --- rewritten 011700

!lll. optimize
      integer lda,n,ipvt(*),info
      real a(lda,*)
!
!     sgefa factors a real matrix by gaussian elimination.
!
!     sgefa is usually called by sgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
!
!     on entry
!
!        a       real(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgesl or sgedi will divide by zero
!                     if called.  use  rcond  in sgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 07/14/77 .
!     cleve moler, university of new mexico, argonne national labs.
!
!     subroutines and functions
!
!     blas saxpy,sscal,isamax
!
!     internal variables
!
      real t
      integer isamax,j,k,kp1,l,nm1

!
! --- gaussian elimination with partial pivoting
!
      info    = 0
      ipvt(n) = n
      if (a(n,n) .eq. 0.0) info = n
      if (n .lt. 2) return

      nm1 = n - 1
      do 60 k = 1, nm1
         kp1 = k + 1

	! --- find l = pivot index
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l

	! --- zero pivot implies this column already triangularized
         if (a(l,k) .eq. 0.0) go to 40

	! --- interchange if necessary
            if (l .ne. k)then
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
	    endif

	! --- compute multipliers
            t = -1.0/a(k,k)
            call sscal(n-k,t,a(k+1,k),1)
	    
	! --- row elimination with column indexing
            do 30 j = kp1, n
               t = a(l,j)
               if (l .ne. k)then
                  a(l,j) = a(k,j)
                  a(k,j) = t
	       endif
               call saxpysmp(n-k,t,a(k+1,k),a(k+1,j))
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue

      return
      end

! ===============================================
! -- TC8
!
      subroutine sgbfa(abd,lda,n,ml,mu,ipvt,info)
!     ===========================================

!lll. optimize
      integer lda,n,ml,mu,ipvt(*),info
      real abd(lda,*)
!
!     sgbfa factors a real band matrix by elimination.
!
!     sgbfa is usually called by sgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!
!     on entry
!
!        abd     real(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgbsl will divide by zero if
!                     called.  use  rcond  in sgbco for a reliable
!                     indication of singularity.
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do 20 j = 1, n
!                      i1 = max0(1, j-mu)
!                      i2 = min0(n, j+ml)
!                      do 10 i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                10    continue
!                20 continue
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     linpack. this version dated 07/14/77 .
!     cleve moler, university of new mexico, argonne national labs.
!
!     subroutines and functions
!
!     blas saxpy,sscal,isamax
!     fortran max0,min0
!
!     internal variables
!
      real t
      integer i,isamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1

      m = ml + mu + 1
      info = 0

! --- zero initial fill-in columns

      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0

! --- gaussian elimination with partial pivoting

      nm1 = n - 1
      if (nm1 .lt. 1) go to 140
      do 130 k = 1, nm1
         kp1 = k + 1

! --- zero next fill-in column

         jz = jz + 1
         if (jz .gt. n) go to 60
            if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0
   40       continue
   50       continue
   60    continue

! --- find l = pivot index

         lm = min0(ml,n-k)
         l = isamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m

! --- zero pivot implies this column already triangularized

         if (abd(l,k) .eq. 0.0) go to 110

	! --- interchange if necessary
            if (l .eq. m) go to 70
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   70       continue

	! --- compute multipliers
            t = -1.0/abd(m,k)
            call sscal(lm,t,abd(m+1,k),1)

	! --- row elimination with column indexing
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 100
            do 90 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 80
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   80          continue
               call saxpysmp(lm,t,abd(m+1,k),abd(mm+1,j))
   90       continue
  100       continue
         go to 120
  110    continue
            info = k
  120    continue
  130 continue
  140 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0) info = n
      return
      end

! ========================================
! -- TC9
!
      subroutine sgesl(a,lda,n,ipvt,b,job)
!     ===================================

! --- Chien Wang
! --- MIT
! --- rewritten 072198

	integer :: lda,n,ipvt(*),job
	real    :: a(lda,*),b(*)
!
!     sgesl solves the real system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by sgeco or sgefa.
!
!     on entry
!
!        a       real(lda, n)
!                the output from sgeco or sgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from sgeco or sgefa.
!
!        b       real(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b  where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if sgeco has set rcond .gt. 0.0
!        or sgefa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call sgesl(a,lda,n,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 07/14/77 .
!     cleve moler, university of new mexico, argonne national labs.
!
!     subroutines and functions
!
!     blas saxpy,sdot
!
!     internal variables
!
	real    :: sdot,t
	integer :: k,kb,l,nm1

	nm1 = n - 1
	if_job: select case (job)
	  
	  ! --- job = 0 , solve  a * x = b
	  case (0)
	    ! --- first solve  l*y = b
	    if (nm1 .ge. 1)then
		do k = 1, nm1
		  l = ipvt(k)
		  t = b(l)
		  if (l .ne. k)then
			b(l) = b(k)
			b(k) = t
		  endif
		  call saxpysmp(n-k,t,a(k+1,k),b(k+1))
		end do
	    endif
	    ! --- now solve  u*x = y
	    do kb = 1, n
		k = n + 1 - kb
		b(k) = b(k)/a(k,k)
		t = -b(k)
		call saxpysmp(k-1,t,a(1,k),b(1))
	    end do

	  ! --- job = nonzero, solve  trans(a) * x = b
	  case default
	    ! --- first solve  trans(u)*y = b
	    do k = 1, n
		t = sdot(k-1,a(1,k),1,b(1),1)
		b(k) = (b(k) - t)/a(k,k)
	    end do
	    ! --- now solve trans(l)*x = y
	    if (nm1 .ge. 1)then
		do kb = 1, nm1
		  k = n - kb
		  b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
		  l = ipvt(k)
		  if (l .ne. k)then
			t = b(l)
			b(l) = b(k)
			b(k) = t
		  endif
		end do
	    endif
	    
	end select if_job

      	return
	 end

! ================================================
! -- TC10
!
      subroutine sgbsl(abd,lda,n,ml,mu,ipvt,b,job)
!     ===========================================

! --- Chien Wang
! --- MIT
! --- rewritten 072198

	integer :: lda,n,ml,mu,ipvt(*),job
	real    :: abd(lda,*),b(*)
!
!     sgbsl solves the real band system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by sgbco or sgbfa.
!
!     on entry
!
!        abd     real(lda, n)
!                the output from sgbco or sgbfa.
!
!        lda     integer
!                the leading dimension of the array  abd .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!
!        mu      integer
!                number of diagonals above the main diagonal.
!
!        ipvt    integer(n)
!                the pivot vector from sgbco or sgbfa.
!
!        b       real(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b , where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if sgbco has set rcond .gt. 0.0
!        or sgbfa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 07/14/77 .
!     cleve moler, university of new mexico, argonne national labs.
!
!     subroutines and functions
!
!     blas saxpy,sdot
!     fortran min0
!
!     internal variables
!
	real    :: sdot,t
	integer :: k,kb,l,la,lb,lm,m,nm1

	m = mu + ml + 1
	nm1 = n - 1
	
	if_job: select case (job)
	  ! --- job = 0 , solve  a * x = b
	  case (0)
	    ! --- first solve l*y = b
	    if (ml .ne. 0 .and. nm1 .ge. 1) then
		do k = 1, nm1
		  lm = min0(ml,n-k)
		  l = ipvt(k)
		  t = b(l)
		  if (l .ne. k)then
			b(l) = b(k)
			b(k) = t
		  endif
		  call saxpysmp(lm,t,abd(m+1,k),b(k+1))
		end do
	    endif
	    ! --- now solve  u*x = y
	    do kb = 1, n
		k = n + 1 - kb
		b(k) = b(k)/abd(m,k)
		lm = min0(k,m) - 1
		la = m - lm
		lb = k - lm
		t = -b(k)
		call saxpysmp(lm,t,abd(la,k),b(lb))
	    end do

	  ! --- job = nonzero, solve  trans(a) * x = b
	  case default
	    ! --- first solve  trans(u)*y = b
	    do k = 1, n
		lm = min0(k,m) - 1
		la = m - lm
		lb = k - lm
		t = sdot(lm,abd(la,k),1,b(lb),1)
		b(k) = (b(k) - t)/abd(m,k)
	    end do
	    ! --- now solve trans(l)*x = y
	    if (ml .ne. 0 .and. nm1 .ge. 1) then
		do kb = 1, nm1
		  k = n - kb
		  lm = min0(ml,n-k)
		  b(k) = b(k) + sdot(lm,abd(m+1,k),1,b(k+1),1)
		  l = ipvt(k)
		  if (l .ne. k) then
			t = b(l)
			b(l) = b(k)
			b(k) = t
		  endif
		end do
	    endif

	end select if_job
	
	return
	 end

! ==================================
! -- TC11
!
      subroutine sscal(n,sa,sx,incx)
!     ==============================

! --- Chien Wang
! --- MIT
! --- rewritten 072198

!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to 1.
!     jack dongarra, linpack, 6/17/77.
!
	real    :: sa,sx(*)
	integer :: i,incx,m,mp1,n,nincx

	if(n.le.0)return
	
	select case (incx)
	  case(1)
	  ! --- code for increment equal to 1
	  ! --- clean-up loop
	  	m = mod(n,5)
	  	if( m .ne. 0 ) then
		  do i = 1,m
		  	sx(i) = sa*sx(i)
		  end do
		  if( n .lt. 5 ) return		
	  	endif
	  	mp1 = m + 1
	  	do i = mp1,n,5
		  sx(i) = sa*sx(i)
		  sx(i + 1) = sa*sx(i + 1)
		  sx(i + 2) = sa*sx(i + 2)
		  sx(i + 3) = sa*sx(i + 3)
		  sx(i + 4) = sa*sx(i + 4)
	  	end do
	  case default
	  ! --- code for increment not equal to 1
	  	nincx = n*incx
	  	do i = 1,nincx,incx
		  sx(i) = sa*sx(i)
	  	end do
	end select

      	return
	 end


! ===================================
! -- TC12
!
      subroutine saxpysmp(n,sa,sx,sy)
!     ===============================

! =====================================================
!	saxpysmp:  a simplified version of saxpy with
!		   incx = incy = 1
!	---------------------------------------------
!	Chien Wang
!	MIT
!
!	Creates:	010795
!	Last revision:	072198
! =====================================================

!
!     constant times a vector plus a vector.
!     uses unrolled loop for increments equal to one.
!     jack dongarra, linpack, 6/17/77.
!
	real    :: sx(*),sy(*),sa
	integer :: i,ix,iy,m,mp1,n

	if(n.le.0)return
	if (sa .eq. 0.0) return

	! --- code for both increments equal to 1
	! --- clean-up loop
	m = mod(n,4)
	if( m .ne. 0 ) then
	  sy(1:m) = sy(1:m) + sa*sx(1:m)
	  if( n .lt. 4 ) return	  
   	endif
	mp1 = m + 1
	do i = mp1,n,4
	  sy(i) = sy(i) + sa*sx(i)
	  sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
	  sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
	  sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
	end do
	
	return
	 end

! ======================================================================
! -- TC13
!
      subroutine xerrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
!     ==================================================================

! --- Chien Wang
! --- MIT
! --- rewritten 072198

      integer ::  nmes, msg(nmes), nerr, level, ni, i1, i2, nr,   &
                 i, lun, lunit, mesflg, ncpw, nch, nwds
      real    :: r1, r2

! ----------------------------------------------------------------------
! subroutines xerrwv, xsetf, and xsetun, as given here, constitute
! a simplified version of the slatec error handling package.
! written by a. c. hindmarsh at llnl.  version of march 30, 1987.
!
! all arguments are input arguments.
!
! msg    = the message (hollerith literal or integer array).
! nmes   = the length of msg (number of characters).
! nerr   = the error number (not used).
! level  = the error level..
!          0 or 1 means recoverable (control returns to caller).
!          2 means fatal (run is aborted--see note below).
! ni     = number of integers (0, 1, or 2) to be printed with message.
! i1,i2  = integers to be printed, depending on ni.
! nr     = number of reals (0, 1, or 2) to be printed with message.
! r1,r2  = reals to be printed, depending on nr.
!
! note..  this routine is machine-dependent and specialized for use
! in limited context, in the following ways..
! 1. the number of hollerith characters stored per word, denoted
!    by ncpw below, is a data-loaded constant.
! 2. the value of nmes is assumed to be at most 60.
!    (multi-line messages are generated by repeated calls.)
! 3. if level = 2, control passes to the statement   stop
!    to abort the run.  this statement may be machine-dependent.
! 4. r1 and r2 are assumed to be in single precision and are printed
!    in e21.13 format.
! 5. the common block /eh0001/ below is data-loaded (a machine-
!    dependent feature) with default values.
!    this block is needed for proper retention of parameters used by
!    this routine which the user can reset by calling xsetf or xsetun.
!    the variables in this block are as follows..
!       mesflg = print control flag..
!                1 means print all messages (the default).
!                0 means no printing.
!       lunit  = logical unit number for messages.
!                the default is 6 (machine-dependent).
! ----------------------------------------------------------------------
! the following are instructions for installing this routine
! in different machine environments.
!
! to change the default output unit, change the data statement below.
!
! for some systems, the data statement below must be replaced
! by a separate block data subprogram.
!
! for a different number of characters per word, change the
! data statement setting ncpw below, and format 10.  alternatives for
! various computers are shown in comment cards.
!
! for a different run-abort command, change the statement following
! statement 100 at the end.
! ----------------------------------------------------------------------

	common /eh0001/ mesflg, lunit

! ----------------------------------------------------------------------
! the following data-loaded value of ncpw is valid for the cdc-6600
! and cdc-7600 computers.
!     data ncpw/10/
! the following is valid for the cray-1 computer.
      data ncpw/8/
! the following is valid for the burroughs 6700 and 7800 computers.
!     data ncpw/6/
! the following is valid for the pdp-10 computer.
!     data ncpw/5/
! the following is valid for the vax computer with 4 bytes per integer,
! and for the ibm-360, ibm-370, ibm-303x, and ibm-43xx computers.
!     data ncpw/4/
! the following is valid for the pdp-11, or vax with 2-byte integers.
!     data ncpw/2/
! ----------------------------------------------------------------------
#if (defined SPMD )
	mesflg = 0
	lunit  = 6
#else
	data mesflg/1/, lunit/6/
#endif

	if (mesflg .ne. 0) then
	! --- get logical unit number.
	  lun = lunit
	! --- get number of words in message.
	  nch = min0(nmes,60)
	  nwds = nch/ncpw
	  if (nch .ne. nwds*ncpw) nwds = nwds + 1
	  ! --- write the message.
	  write (lun, 10) (msg(i),i=1,nwds)
! ----------------------------------------------------------------------
! the following format statement is to have the form
! 10  format(1x,mmann)
! where nn = ncpw and mm is the smallest integer .ge. 60/ncpw.
! the following is valid for ncpw = 10.
! 10  format(1x,6a10)
! the following is valid for ncpw = 8.
  10  format(1x,8a8)
! the following is valid for ncpw = 6.
! 10  format(1x,10a6)
! the following is valid for ncpw = 5.
! 10  format(1x,12a5)
! the following is valid for ncpw = 4.
! 10  format(1x,15a4)
! the following is valid for ncpw = 2.
! 10  format(1x,30a2)
! ----------------------------------------------------------------------
	  if (ni .eq. 1) write (lun, 20) i1
 20   format(6x,"in above message,  i1 =",i10)
	  if (ni .eq. 2) write (lun, 30) i1,i2
 30   format(6x,"in above message,  i1 =",i10,3x,"i2 =",i10)
	  if (nr .eq. 1) write (lun, 40) r1
 40   format(6x,"in above message,  r1 =",e21.13)
	  if (nr .eq. 2) write (lun, 50) r1,r2
 50   format(6x,"in above,  r1 =",e21.13,3x,"r2 =",e21.13)

	endif
	
	if (level .ne. 2) stop
	
	return
	 end

! ===============================
! -- TC14
!
      real function r1mach (idum)
!     ===========================

! --- Chien Wang
! --- MIT
! --- rewritten 072198

	integer :: idum
	real    :: u, comp
	
! ----------------------------------------------------------
! this routine computes the unit roundoff of the machine.
! this is defined as the smallest positive machine number
! u such that  1.0 + u .ne. 1.0
! ----------------------------------------------------------

	u = 1.0
10	u = u*0.5
	comp = 1.0 + u
	if (comp .ne. 1.0) go to 10

	r1mach = u*2.0
	
	return
	 end

! =================================
! -- TC15
!
      real function vnorm (n, v, w)
!     =============================

! --- Chien Wang
! --- MIT
! --- rewritten 072198

! ------------------------------------------------------------------
! this function routine computes the weighted root-mean-square norm
! of the vector of length n contained in the array v, with weights
! contained in the array w of length n..
!   vnorm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )
! ------------------------------------------------------------------
	integer :: n, i
	real    :: v(n), w(n), sum

	sum = 0.0
	do i = 1,n
	  sum = sum + (v(i)*w(i))**2
	end do
	
	vnorm = sqrt(sum/real(n))
	
	return
	 end

! ======================================
! -- TC16
!
      integer function isamax(n,sx,incx)
!     ==================================

! --- Chien Wang
! --- MIT
! --- rewritten 072198

!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 6/17/77.
!
	real    :: sx(*),smax
	integer :: i,incx,ix,n

	isamax = 1
	if(n.le.1) return

	select case (incx)
	  case (1)
	  ! --- code for increment equal to 1
	  	smax = abs(sx(1))
	  	do i = 2,n
		  if(abs(sx(i)).gt.smax) then
		  	isamax = i
		  	smax = abs(sx(i))
		  endif
	  	end do
	
	  case default
	  ! --- code for increment not equal to 1
	  	ix = 1
	  	smax = abs(sx(1))
	  	ix = ix + incx
	  	do i = 2,n
		  if(abs(sx(ix)).gt.smax)then
		  	isamax = i
		  	smax = abs(sx(ix))
		  endif
		  ix = ix + incx
	  	end do
	end select
	
	return
	 end


! =========================================
! -- TC17
!
      real function sdot(n,sx,incx,sy,incy)
!     =====================================

! --- Chien Wang
! --- MIT
! --- rewritten 072198

!
!     forms the dot product of a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 6/17/77.
!
	real    :: sx(*),sy(*),stemp
	integer :: i,incx,incy,ix,iy,m,mp1,n

	stemp = 0.0
	sdot = 0.0
	if(n.le.0)return
	
	if(incx.eq.1.and.incy.eq.1) then
	! --- code for both increments equal to 1
	! --- clean-up loop
	  m = mod(n,5)
	  if( m .ne. 0 ) then
		do i = 1,m
		  stemp = stemp + sx(i)*sy(i)
		end do
	  endif
	  if ( n .ge. 5 ) then
		mp1 = m + 1
		do i = mp1,n,5
		  stemp = stemp + sx(i)*sy(i)          &
     			       + sx(i + 1)*sy(i + 1)   &
     			       + sx(i + 2)*sy(i + 2)   &
     			       + sx(i + 3)*sy(i + 3)   &
     			       + sx(i + 4)*sy(i + 4)
	  	end do
	  endif	    
   	else
	! --- code for unequal increments or equal increments
	! --- 	   not equal to 1
	  ix = 1
	  iy = 1
	  if(incx.lt.0)ix = (-n+1)*incx + 1
	  if(incy.lt.0)iy = (-n+1)*incy + 1
	  do i = 1,n
		stemp = stemp + sx(ix)*sy(iy)
		ix = ix + incx
		iy = iy + incy
	  end do
	endif
	
	sdot = stemp 
	
	return
	 end


! =================================================
! -- TC18
!
      subroutine prepj (neq, y, yh, nyh, ewt, ftem,   & 
                       savf, wm, iwm, f, jac)
!     =============================================

! --- Chien Wang
! --- MIT
! --- 072198

	external f, jac
	integer :: neq(*), nyh, iwm(*), iownd, iowns,               &
     		  icf, ierpj, iersl, jcur, jstart, kflag, 	    &
     		  l, meth, miter, maxord, maxcor, msbp, mxncf,	    &
     		  n, nq, nst, nfe, nje, nqu,			    &
     		  i, i1, i2, ier, ii, j, j1, jj, lenp,		    &
     		  mba, mband, meb1, meband, ml, ml3, mu, np1
       real    :: y(*), yh(nyh,*), ewt(*), ftem(*), savf(*), 	    &
     		  wm(*), rowns, ccmax, el0, h, hmin, hmxi,	    &
     		  hu, rc, tn, uround, con, di, fac, hl0,	    &
     		  r, r0, srur, yi, yj, yjj, vnorm

       common /ls0001/ rowns(209),				    &
     	ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,		    &
     	iownd(14), iowns(6),					    &
     	icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,     &
     	maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
     
! ----------------------------------------------------------------------
! prepj is called by stode to compute and process the matrix
! p = i - h*el(1)*j , where j is an approximation to the jacobian.
! here j is computed by the user-supplied routine jac if
! miter = 1 or 4, or by finite differencing if miter = 2, 3, or 5.
! if miter = 3, a diagonal approximation to j is used.
! j is stored in wm and replaced by p.  if miter .ne. 3, p is then
! subjected to lu decomposition in preparation for later solution
! of linear systems with p as coefficient matrix. this is done
! by sgefa if miter = 1 or 2, and by sgbfa if miter = 4 or 5.
!
! in addition to variables described previously, communication
! with prepj uses the following..
! y     = array containing predicted values on entry.
! ftem  = work array of length n (acor in stode).
! savf  = array containing f evaluated at predicted y.
! wm    = real work space for matrices.  on output it contains the
!         inverse diagonal matrix if miter = 3 and the lu decomposition
!         of p if miter is 1, 2 , 4, or 5.
!         storage of matrix elements starts at wm(3).
!         wm also contains the following matrix-related data..
!         wm(1) = sqrt(uround), used in numerical jacobian increments.
!         wm(2) = h*el0, saved for later use if miter = 3.
! iwm   = integer work space containing pivot information, starting at
!         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band
!         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
! el0   = el(1) (input).
! ierpj = output error flag,  = 0 if no trouble, .gt. 0 if
!         p matrix found to be singular.
! jcur  = output flag = 1 to indicate that the jacobian matrix
!         (or approximation) is now current.
! this routine also uses the common variables el0, h, tn, uround,
! miter, n, nfe, and nje.
! ----------------------------------------------------------------------
	nje = nje + 1
	ierpj = 0
	jcur = 1
	hl0 = h*el0

	miter_flow: select case (miter)
	  case (1)
	  ! --- if miter = 1, call jac and multiply by scalar.
		lenp = n*n
		wm(3:lenp+2) = 0.0
		call jac (neq, tn, y, 0, 0, wm(3), n)
		con = -hl0
		wm(3:lenp+2) = wm(3:lenp+2)*con

	  case (2)
	  ! --- if miter = 2, make n calls to f to approximate j.
		fac = vnorm (n, savf, ewt)
		r0 = 1000.0*abs(h)*uround*real(n)*fac
		if (r0 .eq. 0.0) r0 = 1.0
		srur = wm(1)
		j1 = 2
		do j = 1,n
		  yj = y(j)
		  r = max(srur*abs(yj),r0/ewt(j))
		  y(j) = y(j) + r
		  fac = -hl0/r
		  call f (neq, tn, y, ftem)
		  do i = 1,n
			wm(i+j1) = (ftem(i) - savf(i))*fac
 		  end do
		  y(j) = yj
		  j1 = j1 + n
		end do
		nfe = nfe + n

	  case (3)
		wm(2) = hl0
		r = el0*0.1
		y(1:n) = y(1:n) + r*(h*savf(1:n) - yh(1:n,2))
		call f (neq, tn, y, wm(3))
		nfe = nfe + 1
		do i = 1,n
		  r0 = h*savf(i) - yh(i,2)
		  di = 0.1*r0 - h*(wm(i+2) - savf(i))
		  wm(i+2) = 1.0
		  if (abs(r0) .ge. uround/ewt(i)) then
			if (abs(di) .eq. 0.0) then
			  ierpj = 1
			  exit
			endif
			wm(i+2) = 0.1*r0/di
		  endif
		end do

	  case (4)
	  ! --- if miter = 4, call jac and multiply by scalar.
		ml = iwm(1)
		mu = iwm(2)
		ml3 = ml + 3
		mband = ml + mu + 1
		meband = mband + ml
		lenp = meband*n
		wm(3:lenp+2) = 0.0
		call jac (neq, tn, y, ml, mu, wm(ml3), meband)
		con = -hl0
		wm(3:lenp+2) = wm(3:lenp+2)*con
		
	  case (5)
	  ! --- if miter = 5, make mband calls to f to approximate j.
		ml = iwm(1)
		mu = iwm(2)
		mband = ml + mu + 1
		mba = min0(mband,n)
		meband = mband + ml
		meb1 = meband - 1
		srur = wm(1)
		fac = vnorm (n, savf, ewt)
		r0 = 1000.0*abs(h)*uround*real(n)*fac
		if (r0 .eq. 0.0) r0 = 1.0
		do j = 1,mba
		  do i = j,n,mband
			yi = y(i)
			r = max(srur*abs(yi),r0/ewt(i))
			y(i) = y(i) + r
		  end do
		  call f (neq, tn, y, ftem)
		  do jj = j,n,mband
			y(jj) = yh(jj,1)
			yjj = y(jj)
			r = max(srur*abs(yjj),r0/ewt(jj))
			fac = -hl0/r
			i1 = max0(jj-mu,1)
			i2 = min0(jj+ml,n)
			ii = jj*meb1 - ml + 2
			do i = i1,i2
			  wm(ii+i) = (ftem(i) - savf(i))*fac
			end do	  
		  end do
		end do
		nfe = nfe + mba
		
	end select miter_flow
	
	miter_flow_2: select case (miter)
	  case (1, 2)      
	  	! --- add identity matrix.
		j = 3
		np1 = n + 1
		do i = 1,n
		  wm(j) = wm(j) + 1.0
		  j = j + np1
		end do
		! --- do lu decomposition on p.
		call sgefa (wm(3), n, n, iwm(21), ier)
		if (ier .ne. 0) ierpj = 1
		
	  ! --- if miter = 3, construct a diagonal approximation to j and p.
	  
	  case (4, 5)
	  	! --- add identity matrix.
		ii = mband + 2
		do i = 1,n
		  wm(ii) = wm(ii) + 1.0
		  ii = ii + meband
		end do
		! --- do lu decomposition of p.
		call sgbfa (wm(3), meband, n, ml, mu, iwm(21), ier)
      		if (ier .ne. 0) ierpj = 1
		
	end select miter_flow_2
	
	return
	 end
