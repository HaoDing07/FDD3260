program post_process
!
!
!  Fortran routine for post_processing of LES-MPC results in 2D
!  The routine reads the raw data stored in specified files from
!  OUTPUT and creates a single formatted output file readable
!  by vizualisation softwares
!
!
  integer  :: i, j, k, m, n
  integer  :: ifile, itime, iprof, idiagop, iprof2d, itt, ip
  integer  :: nx, nz, ny, nv, nd, ntot
  integer  :: nyp, nny, js, je
  integer  :: iform, skip(100)
  integer  :: ierr=1,lfirst
  integer, allocatable, dimension(:) :: ind, idiag
    
  double precision, allocatable, dimension(:,:,:,:) :: a, d
  double precision, allocatable, dimension(:,:)     :: a1d   
  double precision, allocatable, dimension(:,:,:)   :: u, v, w  
  double precision  :: tt, stime
  real, allocatable, dimension(:,:,:) :: xx, yy, zzz
  real, allocatable, dimension(:)     :: zz
  real     :: dx, dz, dy

  integer  :: in, iq, iqc, iz, ndiam, npts, nran, kk
  real  :: b, am, av, betam, betav
  real  :: ddiam, diam, zm
  real, allocatable, dimension(:) :: diam2, fdm
  real, allocatable, dimension(:,:,:) :: fd
  real, allocatable, dimension(:,:) :: fd2
  integer    :: SEED
  integer    :: tseed(8)
  integer, allocatable, dimension(:)   :: ranx, rany
  
  character(len=40)  :: filename
  character(len=11)  :: form
  character(len=7)   :: dirname='OUTPUT/'
  character(len=13)  :: car4
  character(len=1)   :: car1
  character(len=2)   :: car2
  character(len=3)   :: alt, icename
  character(len=13), allocatable, dimension(:) :: namev
  character(len=6),  allocatable, dimension(:) :: time
!  
  logical  :: l, ex
!
!  Dimensions
!
  print*, 'number of points in x and z:'
  read*, nx, ny, nz
  if (ny <= 4) ny = 1
  print*, 'cell size in x and z:'
  read*, dx, dy, dz
  print*, 'format: 0-binary, 1-ASCII'
  read*, iform
  print*, 'ice species (DEN, ASS, PLA, COL, BUL, ..., NO):'
  read*, icename
  print*, 'altitude slice for size distributions:'
  read*, zm
!  
  if (iform == 1) then
    form = 'formatted'
  else
    form = 'unformatted'
  endif
!
  allocate(xx(1:nx,1:ny,1:nz), zz(1:nz), yy(1:nx,1:ny,1:nz), zzz(1:nx,1:ny,1:nz))
!
!  Default grid dimensions: zz may be overwritten if file Z0 is present
!
  xx(1,:,:) = 0.0
  zz(1)     = 0.0
  yy(:,1,:) = 0.0
  do i = 2, nx
    xx(i,:,:) = xx(i-1,:,:) + dx
  enddo
  do i = 2, nz
    zz(i) = zz(i-1) + dz
  enddo 
  if (ny /= 1 ) then 
    do i = 2, ny
      yy(:,i,:) = yy(:,i-1,:) + dy
    enddo  
  endif
!
!  Files to read
!
  open(10,FILE='file_list.inp')
  
  nv=0
  ntot=0
  itime = 1
  skip = 0
  do while(ierr >= 0)
    ntot = ntot+1
    read(10,*,iostat=ierr) car4
    if (car4(1:1) /= '#') then
      nv = nv+1
    else
      skip(ntot-1) = 1
    endif
!    
    if (car4 == 'time') then
      do while (trim(car4) /= 'name')
        read(10,*,iostat=ierr) car4
	if (trim(car4) /= '0000') itime = itime + 1
      enddo 
      itime = itime - 2
    endif
  enddo
  nv = nv-2
  ntot = ntot-2
  allocate(namev(1:nv))
  
  rewind(10)
  
  if (itime >= 1) then
    allocate(time(1:itime))
    read(10,*) car4
    read(10,*) car4
    if (trim(car4) == '0000') read(10,*) car4
    do i = 1, itime  
      time(i) = car4
      read(10,*) car4
    enddo
  endif
  
  j = 0
  do i = 1, ntot
    if(skip(i) /= 1) then
      j = j+1
      read(10,*) namev(j)
    else
      read(10,*) car4
    endif
  enddo
 
  close(10)
!
!  Open raw data files in OUTPUT, loop over itime
!
  k = 1
  do j = 1, itime
!
  allocate(a(1:nx,1:ny,1:nz,1:nv))
  allocate(a1d(1:nz,1:nv))
  allocate(ind(1:nv), idiag(1:nv))
!
  if (ny == 1) then
    nyp = 1
    nny = 1
    js = 1
    je = 1
  else
    nyp = ny+1
    nny = ny-3
    js = 2
    je = ny-2
  endif
!
!  Read all other files
!
  iprof = 0
  ind = 0
  print*, ''
  do i = 1, nv
    ex  = .true.
    ifile = 11+i
    ip = index(namev(i),'pro')
!
    if (ip == 0) then
!
    if (trim(namev(i)) /= 'Z0') then
      filename = dirname//trim(namev(i))//'.'//trim(time(j))
    else
      filename = dirname//trim(namev(i))
    endif
!
    inquire (FILE=trim(filename), EXIST=ex)
    if (.not. ex) then
      print*, 'File ', trim(filename), ' does not exist'
    endif
!
    if (ex) then
!
    if (filename(7:10) /= '0000') then
    print*, 'Open and read in file ', trim(filename)
    open(ifile,FILE=trim(filename),FORM=form,STATUS='OLD')
!
    if (j == 1 .and. trim(namev(i)) == 'Z0') then
      if (iform == 1) then
        read(ifile,100) zz(:)
      else
        read(ifile) zz(:)      
      endif
      print*, '    - Max height:', maxval(zz(:))  
    else if (ip /= 0) then      
      iprof = iprof+1
      ind(iprof) = i
      if (iform == 1) then
        read(ifile,100) a1d(:,iprof)      
      else
        read(ifile) a1d(:,iprof)    
      endif
      print*, '    - Vertical profile'    
    else if (trim(namev(i)) == 'U') then
      allocate(u(1:nx+1,1:nyp,1:nz+1))
      if (iform == 1) then
        read(ifile,100) u(:,:,:)
      else
      read(ifile) u(:,:,:)      
      endif      
      print*, '    - Max reported value in the file:', maxval(u(:,:,:)), maxloc(u(:,:,:))     
    else if (trim(namev(i)) == 'V') then
      allocate(v(1:nx+1,1:nyp,1:nz+1))
      if (iform == 1) then
        read(ifile,100) v(:,:,:)        
      else
        read(ifile) v(:,:,:)      
      endif
      print*, '    - Max reported value in the file:', maxval(v(:,:,:)), maxloc(v(:,:,:))
    else if (trim(namev(i)) == 'W') then
      allocate(w(1:nx+1,1:nyp,1:nz+1))
      if (iform == 1) then
        read(ifile,100) w(:,:,:)
      else      
        read(ifile) w(:,:,:)
      endif      
      print*, '    - Max reported value in the file:', maxval(w(:,:,:)), maxloc(w(:,:,:))   
    else     
      if (iform == 1) then
        read(ifile,100) a(:,:,:,i)
      else
        read(ifile) a(:,:,:,i)      
      endif
      print*, '    - Max reported value in the file:', maxval(a(:,:,:,i)), maxloc(a(:,:,:,i))
    endif
    endif
!
    close(ifile) 
    endif 
    endif
!
    enddo
!
    do k = 1, nx
      do i = 1, ny
        zzz(k,i,:) = zz(:)
      enddo
    enddo
!
!
!  Calculate derived variables from 2D/3D field
!
    if (j == itime) then
!
!  Define ice microphysics
!
      nd = 3
      allocate(d(1:nx,1:ny,1:nz,1:nd))
!
      b = 2
      if (adjustl(adjustr(icename)) == 'DEN') then
        am = 0.023
	av = 5.02
	betam = 2.29
	betav = 0.48
      else if (adjustl(adjustr(icename)) == 'PLA') then
        am = 1.43
	av = 17.9
	betam = 2.79
	betav = 0.62
      else if (adjustl(adjustr(icename)) == 'ASS') then
        am = 0.0444
	av = 11.72
	betam = 2.1
	betav = 0.41
      else if (adjustl(adjustr(icename)) == 'COL') then	
        am = 0.0093
	av = 30.1
	betam = 1.8
	betav = 0.488
      else if (adjustl(adjustr(icename)) == 'BUL') then	
        am = 0.0533
	av = 220.
	betam = 2.27
	betav = 0.8
      endif
!  	
      do i = 1, nv
        if (trim(namev(i)) == 'QI') then
	  iq = i
	else if (trim(namev(i)) == 'NI') then
	  in = i
	else if (trim(namev(i)) == 'QC') then
	  iqc = i
	endif
      enddo
!
!  First order ice microphysics
!
      d(:,:,:,1) = 0.0
      where(a(:,:,:,iq) > 1.e-10)		&
          d(:,:,:,1) = (cal_gamma(b+betam+1) * am * a(:,:,:,in) / (cal_gamma(b+1)*a(:,:,:,iq)))**(1./betam)
      print*, '    - Computed derived variable: Ice crystals Lambda'
!	  
      d(:,:,:,2) = 0.0
      where(d(:,:,:,1) > 1.e-10)		&
          d(:,:,:,2) = (b+1)/d(:,:,:,1)
      print*, '    - Computed derived variable: Ice crystals mean diameter'
!	  
      d(:,:,:,3) = 0.0
      where(a(:,:,:,iq) > 1.e-10 .and. a(:,:,:,iqc) > 1.e-10)		&
          d(:,:,:,3) = 10.*log10((4.5755e-6*a(:,:,:,in)*(a(:,:,:,iq)/a(:,:,:,in))**2. + 	&
	      1.0785e-6*200.e6*(a(:,:,:,iqc)/200.e+6)**2.) / 2.9568e-18)
      print*, '    - Computed derived variable: Ice crystals mean diameter'
!
!  Size distributions
!
      nran = 20
      ndiam = 3000.
      ddiam = 5e-6
      do k = 2, nz
        if (zz(k) > zm .and. zz(k-1) <= zm) iz = k
      enddo
      allocate(diam2(1:ndiam), fdm(1:ndiam), fd2(1:ndiam,1:nran), fd(1:ndiam,1:nx,1:ny))
      allocate(ranx(1:nran), rany(1:nran))
!
      print*, '    - Compute ice size distributions at z =', floor(zm)
      do n = 1, ny
        do m = 1, nx
          diam = 1.e-5
          do k = 1, ndiam
	    if (m==1 .and. n==1) diam2(k) = diam
	    fd(k,m,n) = 1.1513 * a(m,n,iz,in) * d(m,n,iz,1)**(b+1) * diam**(b+1) * exp(-diam*d(m,n,iz,1))
	    diam = diam + ddiam
	  enddo
	enddo
      enddo
!
      CALL DATE_AND_TIME(VALUES = tseed)
      SEED = tseed(1)+70*(tseed(2)+12*(tseed(3)+31*(tseed(5)+23*(tseed(6)+59*tseed(7)))))
      IF (MOD(SEED,2).EQ.0) SEED = SEED-1
      do k = 1, nran
        ranx(k) = floor(real(nx)*RANF (SEED))
        rany(k) = floor(real(ny)*RANF (SEED))
      enddo
!
      npts = 0
      kk   = 1
      fdm  = 0.
      fd2  = 0.
      do n = 1, ny
        do m = 1, nx
	  npts = npts+1
	  fdm(1:ndiam) = fdm(1:ndiam) + fd(1:ndiam,m,n)
	  do k = 1, nran
	    if (m == ranx(k) .and. n == rany(k) .and. kk <= nran) then
	      fd2(:,kk) = fd(:,m,n)
	      kk = kk + 1
	    endif
	  enddo
        enddo
      enddo
      fdm = fdm/real(npts)
!
!  Write derived variables
!
      write(alt,'(i3)') floor(zm)
      filename = 'size_distribution.'//trim(time(j))//'.'//trim(adjustl(alt))//'.dat'
      open(12,FILE=trim(filename),FORM='formatted')
!
      do k = 1, ndiam
        write(12,*) 1.e+6*diam2(k), fdm(k)/1.e+3, fd2(k,:)/1.e+3
      enddo
    endif
!
100  format(6e15.7)
!
!  Terminate
!  
  if(allocated(a)) deallocate(a)
  if(allocated(a1d)) deallocate(a1d)
  if(allocated(ind)) deallocate(ind)
  if(allocated(idiag)) deallocate(idiag)
  if(allocated(u)) deallocate(u)
  if(allocated(v)) deallocate(v)
  if(allocated(w)) deallocate(w)
  inquire(UNIT=ifile, OPENED=l)
  if (l) close(ifile)
  inquire(UNIT=10, OPENED=l)
  if (l) close(10)
  inquire(UNIT=11, OPENED=l)
  if (l) close(11)
  inquire(UNIT=12, OPENED=l)
  if (l) close(12)
!
  enddo
!
end program

! ============================
	
  function cal_gamma (x)

  IMPLICIT NONE

  real     :: x, xx, cal_gamma
  real     :: pi=3.14159
  real, parameter :: p0=1.000000000190015,	  &
  		     p1=76.18009172947146,	  &
  		     p2=-86.50532032941677,	  &
  		     p3=24.01409824083091,	  &
  		     p4=-1.231739572450155,	  &
  		     p5=1.208650973866179e-3,	  &
  		     p6=-5.395239384953e-6

  xx = p0 + p1/(x + 1.) + p2/(x + 2.) + p3/(x + 3.) + p4/(x + 4.) + p5/(x + 5.) + p6/(x + 6.)
  cal_gamma = sqrt(2.*pi) * xx/x * (x + 5.5)**(x+0.5) * exp(-(x + 5.5))

  return
  end function    


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
