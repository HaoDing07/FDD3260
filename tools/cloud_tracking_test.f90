
  program cloud_tracking
!      =================================================                            

IMPLICIT NONE

  integer :: nx=64, ny=64, nens=20
  real :: dx=500., dy=500.
!
  integer :: i, j, m, n, l, iens, nclouds, icloud, maxs, maxs_t
  real :: nav, nvar, nav_t, nvar_t, nclouds_t, di, dj
  character(len=40) :: filename
  logical :: ex,merge,done
!    
  INTEGER, DIMENSION (8) :: T
  REAL :: RANF,tmp,pi=3.14159265
  INTEGER :: SEED      
!     
  integer, allocatable, dimension(:,:) :: ifdrop, id_cl
  real, allocatable, dimension(:,:) :: ql, xx, yy, dist
  real, allocatable, dimension(:) :: pdf
!
  integer, dimension(10000) :: npts
  real, dimension(10000) :: ic, jc, dav, dvar
  integer, dimension(10000,100) :: ii, jj
  real, dimension(10000,10000) :: dd
  real, dimension(100) :: pdf_t, pdf_d, pdf_dt
!
!----------------------------------------------------------!
!                    Initializations                       !
!----------------------------------------------------------!
!
  allocate(ifdrop(1:nx,1:ny), id_cl(1:nx,1:ny))
  allocate(ql(1:nx,1:ny), xx(1:nx,1:ny), yy(1:nx,1:ny), dist(1:nx,1:ny))
!
!  Ensemble
!
  nav_t=0.
  nvar_t=0.
  pdf_t=0.
  pdf_dt=0.
  nclouds_t=0
  maxs_t=0.
  do iens = 1, nens      
!
  nclouds=0
  id_cl = 0
  ifdrop = 0
  dist = 0.
  ql = 0.
  dav = 0.
  dvar = 0.
!
  do n = 1, 10000
    npts(n) = 0
    dd(n,:) = 0.
    ii(n,:) = 0
    jj(n,:) = 0
    ic(n) = 0
    jc(n) = 0
  enddo
!
!  Random clouds
!
  CALL DATE_AND_TIME(VALUES = T)
  SEED = T(1)+(70+iens)*(T(2)+12*(T(3)+31*(T(5)+23*(T(6)+59*T(7)))))
  IF (MOD(SEED,2).EQ.0) SEED = SEED-1
!
  do j=1,ny
    do i=1,nx
	  tmp = RANF(SEED)
	  xx(i,j) = dx*(i-1)
	  yy(i,j) = dy*(j-1)
!
	  ql(i,j) = max(0.,5.e-3*(tmp-0.8))
      if (ql(i,j) > 1.e-6) ifdrop(i,j) = 1
    enddo
  enddo
!
!----------------------------------------------------------!
!                  Starting cloud search                   !
!----------------------------------------------------------!
!
      do i = 1, nx
        do j = 1, ny
          if ( ifdrop(i,j) == 1 ) then
!
!  Cloudy grid cell detected: does it belong to an existing cloud?
!
			ex=.false.
			if ( ifdrop(i,max(j-1,1)) == 1 .and. j /= 1 ) then
              do n = 1, nclouds
                do m = 1, npts(n)
			      if (i == ii(n,m) .and. j-1 == jj(n,m)) then
			        npts(n)=npts(n)+1
			        icloud=n
			        ex=.true.
			        exit
			      endif
			    enddo
                if (ex) exit
              enddo
			else if ( ifdrop(max(i-1,1),j) == 1 .and. i /= 1 ) then
              do n = 1, nclouds
                do m = 1, npts(n)
			      if (i-1 == ii(n,m) .and. j == jj(n,m)) then
			        npts(n)=npts(n)+1
			        icloud=n
			        ex=.true.
			        exit 
			      endif
                enddo
                if (ex) exit
              enddo
!			else if ( ifdrop(max(i-1,1),max(j-1,1)) == 1 .and. (i /= 1 .or. j /= 1) ) then
!              do n = 1, nclouds
!                do m = 1, npts(n)
!                  if (i-1 == ii(n,m) .and. j-1 == jj(n,m)) then
!			        npts(n)=npts(n)+1
!			        icloud=n
!			        ex=.true.
!			        exit
!			      endif
!                enddo
!              enddo
			endif    
!
			if (.not.ex) then
			  nclouds=nclouds+1
			  npts(nclouds)=1
			  icloud=nclouds
			endif
!
			ii(icloud,npts(icloud))=i
			jj(icloud,npts(icloud))=j
			id_cl(i,j)=icloud
          endif
        enddo
      enddo
!
!----------------------------------------------------------!
!                  Merge detected clouds                   !
!----------------------------------------------------------!
!
 	  do n = 1, nclouds
 	    do m = 1, nclouds
!
 	      merge=.false.
 	      if (m /= n) then
   	        do l = 1, npts(n)
 	          if ( (any(ii(n,l)-1 == ii(m,:) .and. jj(n,l) == jj(m,:))) .or.					&
 	          	   (any(jj(n,l)-1 == jj(m,:) .and. ii(n,l) == ii(m,:))) .or.					&
 	          	   (ii(n,l) == 1 .and. any(ii(m,:) == nx) .and. any(jj(n,l) == jj(m,:)))  		&
! 	          	   .or. jj(n,l)-1 == jj(m,:) .or. jj(n,l)+1 == jj(m,:))) 					&
 	          	   .or. (jj(n,l) == 1 .and. any(jj(m,:) == ny) .and. any(ii(n,l) == ii(m,:))) ) then	
! 	          	   .or. ii(n,l)-1 == ii(m,:) .or. ii(n,l)+1 == ii(m,:))) ) then
 	            merge=.true.
 	            exit
 	          endif
 	        enddo
!
!  Merging 2 clouds
!
 	        if (merge) then
   	          do l = 1, npts(n)
 	            id_cl(ii(n,l),jj(n,l)) = id_cl(ii(m,1),jj(m,1))
   	          enddo
!
   	          ii(n,npts(n)+1:npts(n)+npts(m))=ii(m,1:npts(m))
 	          jj(n,npts(n)+1:npts(n)+npts(m))=jj(m,1:npts(m))
 	          npts(n)=npts(n)+npts(m)
!
		  ii(m,:) = 0
		  jj(m,:) = 0
		  npts(m) = 0
 	        endif
 	      endif
!
 	    enddo
 	  enddo
!
!  Recount clouds
!
  m=0
  do i = 1, nclouds
    if (any(id_cl == i)) then
      m=m+1
    endif
  enddo
  nclouds_t=nclouds_t+m
!
!----------------------------------------------------------!
!                    Find cloud centres     	  	   !
!----------------------------------------------------------!
!
  do n = 1, nclouds
    if (.not.((any(ii(n,:) == 1) .and. any(ii(n,:) == nx)) .or. (any(jj(n,:) == 1) .and. any(jj(n,:) == ny)))) then
      ic(n) = sum(real(ii(n,:)))
      jc(n) = sum(real(jj(n,:)))
    else if ((any(jj(n,:) == 1) .and. any(jj(n,:) == ny)) .and. (any(ii(n,:) == 1) .and. any(ii(n,:) == nx))) then
      do i = 1, npts(n)
        if (jj(n,i) <= ny/2) then
          jc(n) = jc(n) + real(jj(n,i))
        else
          jc(n) = jc(n) + real(jj(n,i)-ny)
        endif
        if (ii(n,i) <= nx/2) then
          ic(n) = ic(n) + real(ii(n,i))
        else
          ic(n) = ic(n) + real(ii(n,i)-nx)
        endif
      enddo
    else if (any(jj(n,:) == 1) .and. any(jj(n,:) == ny)) then
      do i = 1, npts(n)
        if (jj(n,i) <= ny/2) then
          jc(n) = jc(n) + real(jj(n,i))
        else
          jc(n) = jc(n) + real(jj(n,i)-ny)
        endif
        ic(n) = ic(n) + real(ii(n,i))
      enddo
    else if (any(ii(n,:) == 1) .and. any(ii(n,:) == nx)) then
      do i = 1, npts(n)
        if (ii(n,i) <= nx/2) then
          ic(n) = ic(n) + real(ii(n,i))
        else
          ic(n) = ic(n) + real(ii(n,i)-nx)
        endif
        jc(n) = jc(n) + real(jj(n,i))
      enddo
    endif
  enddo
  ic = ic / real(npts)
  jc = jc / real(npts)
!
!----------------------------------------------------------!
!		   	  Calculate distance between clouds			   !
!----------------------------------------------------------!
!
  pdf_d=0.0
  done=.true.
  do n = 1, nclouds
    if (npts(n) > 0) then
      do i = 1, nclouds
        if (npts(i) > 0 .and. i /= n) then
          di = min(abs(ic(i) - ic(n)),abs(abs(ic(i) - ic(n))-nx))
          dj = min(abs(jc(i) - jc(n)),abs(abs(jc(i) - jc(n))-ny))
          dd(n,i) = sqrt((dx*real(di))**2. + (dy*real(dj))**2.)
          do l = 1, 100
            if (dd(n,i) > real(l-1)*dx*nx/50. .and. dd(n,i) <= real(l)*dx*nx/50. .and. dd(n,i) < dx*nx/2.) then
              pdf_d(l) = pdf_d(l)+1
              exit
            endif
          enddo
!
          if (done) then
            do j = 1, npts(i)
              dist(ii(i,j),jj(i,j)) = dd(n,i)
            enddo
          endif
        endif
      enddo
      done=.false.
!
      l=0
      do j = 1, nclouds
        if (dd(n,j) < dx*nx/2.) then
          dav(n) = dav(n) + dd(n,j)
          dvar(n) = dvar(n) + dd(n,j)*dd(n,j)
          l=l+1
        endif
      enddo
      dav(n) = dav(n)/real(l)
      dvar(n) = sqrt(dvar(n)/real(l) - dav(n)**2.)
    endif
  enddo
  pdf_dt=pdf_dt+pdf_d
!
!----------------------------------------------------------!
!                  		Cloud stats     	  	           !
!----------------------------------------------------------!
!
  nav  = sum(real(npts))/real(m)
  nvar = sqrt(sum(real(npts*npts))/real(m) - nav**2.)
  maxs = maxval(npts)
  maxs_t=max(maxs_t,maxs)
  print*, iens, m, nav, nvar
!
  allocate(pdf(1:maxs))
  pdf=0.0
  do n = 1, nclouds
    do m = 1, maxs
      if (npts(n) == m) then
        pdf(m)=pdf(m)+1
        exit
      endif
    enddo
  enddo
!
!  Ensemble statistiis
!
  nav_t=nav_t+nav 
  nvar_t=nvar_t+nvar
  pdf_t(1:maxs)=pdf_t(1:maxs)+pdf
  deallocate(pdf)
!
  enddo
!
  nclouds_t=real(nclouds_t)/real(nens)
  nav_t=nav_t/real(nens)
  nvar_t=nvar_t/real(nens)
  pdf_t=pdf_t/real(nens)
  pdf_dt=pdf_dt/real(nens)
!
  open(11,FILE='cloud_stats.csv')
  write(11,*) 'size,    pdf,    pdf_d'
  do i = 1, 25
    write(11,*) real(i), ',', pdf_t(i), ',', pdf_dt(i)
  enddo
  close(11)
!
!----------------------------------------------------------!
!                  			Output     	  		           !
!----------------------------------------------------------!
!
  filename = 'cloud_tracking_test.tec'
  open(10,FILE=trim(filename),FORM='formatted')
!
  write(10,*) 'TITLE     = "output"'
  write(10,*) 'VARIABLES = "x"'
  write(10,*) '"y"'
  write(10,*) '"Ql"'
  write(10,*) '"Mask"'
  write(10,*) '"ID"'
  write(10,*) '"Dist"'
!
  write(10,*) 'ZONE T="Zone"'
  write(10,*) 'I=',nx, ', J=', ny, ', K=', 1
  write(10,*) 'DATAPACKING=BLOCK'
!
  write(10,101) xx(1:nx,1:ny)
  write(10,101) yy(1:nx,1:ny)
  write(10,101) ql(1:nx,1:ny)
  write(10,101) real(ifdrop(1:nx,1:ny))
  write(10,101) real(id_cl(1:nx,1:ny))
  write(10,101) dist(1:nx,1:ny)
!
101  format(3(e15.7,3x))  
!
  close(10)
!
  deallocate (ifdrop, id_cl, ql, xx, yy, dist)
!
!----------------------------------------------------------!
!
  return
  end program cloud_tracking


FUNCTION RANF(SEED) RESULT (CR)
!
! Function to generate a uniform random number in [0,1]
! following x(i+1)=a*x(i) mod c with a=7** 5 and
! c=2** 31-1.  Here the seed is a global variable.
!
  IMPLICIT NONE
  INTEGER :: H, L, T, A, C, Q, R, SEED
  DATA A/16807/, C/2147483647/, Q/127773/, R/2836/
  REAL :: CR
!
  H = SEED/Q
  L = MOD(SEED, Q)
  T = A*L - R*H
  IF (T .GT. 0) THEN
    SEED = T
  ELSE
    SEED = C + T
  END IF
  CR = SEED/FLOAT(C)
!  CR = 1. - 2.*CR
END FUNCTION RANF
