program post_process
!
!
!  Fortran routine for post_processing of LES-MPC results in 2D
!  The routine reads the raw data stored in specified files from
!  OUTPUT and creates a single formatted output file readable
!  by vizualisation softwares
!
!
  integer  :: i, j, k, l, m, n, nn, mm
  integer  :: ifile, itime, iprof, idiagop, iprof2d, itt, ispro, ip, ideriv, iseries, i3d
  integer  :: nx, nz, ny, nv, ntot, sumid
  integer  :: nyp, nny, js, je
  integer  :: iform, skip(100)
  integer  :: ierr=1,lfirst
  integer  :: ivv, inr, iqr, ini, iqi, iqid, inc, iqc, indf, inbcf, irad, item
  integer, allocatable, dimension(:) :: ind, idiag
    
  double precision, allocatable, dimension(:,:,:,:) :: a  
  double precision, allocatable, dimension(:,:)     :: a1d   
  double precision, allocatable, dimension(:,:,:)   :: a2d   
  double precision, allocatable, dimension(:,:,:,:) :: a2dd, ad1d2
  double precision, allocatable, dimension(:,:)     :: d1d2, series
  double precision, allocatable, dimension(:)       :: d1du, d1dv, d1dw, d1d1, vprecr, vpreci
  double precision, allocatable, dimension(:,:,:)   :: u, v, w  
  double precision  :: tt, stime, sum
  real, allocatable, dimension(:)     :: zz
  real, allocatable, dimension(:,:,:) :: xx, yy, zzz
  real     :: dx, dz, dy
  
  character(len=40)  :: filename
  character(len=11)  :: form
  character(len=7)   :: dirname='OUTPUT/'
  character(len=13)  :: car4
  character(len=1)   :: car1
  character(len=2)   :: car2
  character(len=4)   :: car3
  character(len=13), allocatable, dimension(:) :: namev, namev2
  character(len=6),  allocatable, dimension(:) :: time, time2
!  
  logical  :: ex, ex2, ex4, ex5, ll
!
!  Dimensions
!
  print*, 'number of points in x and z:'
  read*, nx, ny, nz
  print*, 'cell size in x and z:'
  read*, dx, dy, dz
  print*, 'Profiles format: 0-1d or 1-2d'
  read*, iprof2d
  print*, 'Output 3d fields?'
  read*, i3d
  print*, 'Read diagnostics?'
  read*, idiagop
  print*, 'Calculate derived variables?'
  read*, ideriv
  print*, 'Calculate derived time-series?'
  read*, iseries
  print*, 'Start time for diagnostics and 2d profiles'
  read*, stime
  print*, 'format: 0-binary, 1-ASCII'
  read*, iform
!  
  if (iform == 1) then
    form = 'formatted'
  else
    form = 'unformatted'
  endif
!
  if (ny <= 4) then
    ny = 1
  endif
!
  allocate(xx(1:nx,1:ny,1:nz), yy(1:nx,1:ny,1:nz), zz(1:nz), zzz(1:nx,1:ny,1:nz))
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
    nyp = ny
    nny = ny
    js = 1
    je = ny
  endif
!
!  Read all other files
!
  ini = 0; inr = 0; iqi = 0; iqr = 0; ivv = 0
  inbcf = 0.; indf = 0.
  iprof = 0
  ind = 0
  ex2 = .true.
  print*, ''
  do i = 1, nv
    ex  = .true.
    ifile = 11+i
    ip = index(namev(i),'pro')

    if (trim(namev(i)) /= 'Z0') then
      filename = dirname//trim(namev(i))//'.'//trim(time(j))
    else
      filename = dirname//trim(namev(i))
    endif
!
    inquire (FILE=trim(filename), EXIST=ex)
    if (.not. ex) then
      print*, 'File ', trim(filename), ' does not exist'
      if (filename(1:1) == 'U') then
        if (filename(3:5) == 'pro') ex2 = .false.
      endif
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
      print*, '    - Max height:', maxval(zz)  
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
      allocate(u(1:nx,1:nyp,1:nz))
      if (iform == 1) then
        read(ifile,100) u(:,:,:)
      else
      read(ifile) u(:,:,:)      
      endif      
      print*, '    - Max reported value in the file:', maxval(u(:,:,:)), maxloc(u(:,:,:))     
    else if (trim(namev(i)) == 'V' .and. ny /= 1) then
      ivv = 1
      allocate(v(1:nx,1:nyp,1:nz))
      if (iform == 1) then
        read(ifile,100) v(:,:,:)        
      else
        read(ifile) v(:,:,:)      
      endif
      print*, '    - Max reported value in the file:', maxval(v(:,:,:)), maxloc(v(:,:,:))
    else if (trim(namev(i)) == 'W') then
      allocate(w(1:nx,1:nyp,1:nz))
      if (iform == 1) then
        read(ifile,100) w(:,:,:)
      else      
        read(ifile) w(:,:,:)
      endif      
      print*, '    - Max reported value in the file:', maxval(w(:,:,:)), maxloc(w(:,:,:))   
    else if (trim(namev(i)) /= 'Z0') then 
      if (iform == 1) then
        read(ifile,100) a(:,:,:,i)
      else
        read(ifile) a(:,:,:,i)      
      endif
      print*, '    - Max reported value in the file:', maxval(a(:,:,:,i)), maxloc(a(:,:,:,i))
      if (trim(namev(i)) == 'QC') iqc = i
      if (trim(namev(i)) == 'NC') inc = i
      if (trim(namev(i)) == 'QR') iqr = i
      if (trim(namev(i)) == 'NR') inr = i
      if (trim(namev(i)) == 'QI') iqi = i
      if (trim(namev(i)) == 'NI') ini = i
      if (trim(namev(i)) == 'INDF') indf = i
      if (trim(namev(i)) == 'INBCF') inbcf = i
      if (trim(namev(i)) == 'RADDT') irad = i
    endif
    endif
!
    close(ifile) 
    endif 
!
  enddo
!
  do k = 1, nx
    do l = 1, ny
      zzz(k,l,:) = zz(:)
    enddo
  enddo
!
!  Write 2D/3D data in output
!
    if (ex2 .and. i3d == 1) then
      filename = 'les_output.'//trim(time(j))//'.tec'
      open(10,FILE=trim(filename),FORM='formatted')
!
      write(10,*) 'TITLE     = "output"'
      write(10,*) 'VARIABLES = "x"'
      write(10,*) '"y"'
      write(10,*) '"z"'
      do i = 1, nv
        if (trim(namev(i)) /= 'Z0' .and. trim(namev(i)) /= 'K' .and. index(namev(i),'pro') == 0)    &
          write(10,*) '"'//trim(namev(i))//'"'
        if (trim(namev(i)) == 'K')     &
          write(10,*) '"Ksgs"'
      enddo
      write(10,*) 'ZONE T="Zone"'
      write(10,*) 'I=',nx, ', J=', ny, ', K=', nz-1
      write(10,*) 'DATAPACKING=BLOCK'
!
      write(10,101) xx(1:nx,js:je,1:nz-1)
      write(10,101) yy(1:nx,js:je,1:nz-1)
      write(10,101) zzz(1:nx,js:je,1:nz-1)
      do i = 1, nv
        if (trim(namev(i)) == 'U') then
        write(10,101) real(u(1:nx,js:je,1:nz-1),4)
        else if (trim(namev(i)) == 'V' .and. ny /= 1) then
        write(10,101) real(v(1:nx,js:je,1:nz-1),4)
        else if (trim(namev(i)) == 'W') then
        write(10,101) real(w(1:nx,js:je,1:nz-1),4)                
        else if (trim(namev(i)) /= 'Z0' .and. index(namev(i),'pro') == 0)then
        write(10,101) real(a(1:nx,js:je,1:nz-1,i),4)
        endif
      enddo 
    endif
!
!  Calculate and write derived variables in 2D (z - t)
!
    if (ex2 .and. ideriv == 1) then        
      if (j == 1) allocate(d1d1(1:nz), d1d2(1:nz,1:105), d1du(1:nz), d1dv(1:nz), d1dw(1:nz), 		&
          namev2(1:15), vprecr(1:nz), vpreci(1:nz), ad1d2(1:nx,1:ny,1:nz,1:45))	
!
      d1d2 = 0.
      mm = 0
      do i = 1, nv      
        if (trim(namev(i)) == 'U') then
	  nn = 0
	  d1du = 0.0
      do l = js, je
	    do k = 1, nx
	      nn = nn + 1
	      d1du(1:nz) = d1du(1:nz) + u(k,l,1:nz)
	    enddo
	  enddo 
	  d1du = d1du / real(nn)     
        else if (trim(namev(i)) == 'V' .and. ny /= 1) then
	  nn = 0
	  d1dv = 0.0
          do l = js, je
	    do k = 1, nx
	      nn = nn + 1
	      d1dv(1:nz) = d1dv(1:nz) + v(k,l,1:nz)
	    enddo
	  enddo 
	  d1dv = d1dv / real(nn)    
        else if (trim(namev(i)) == 'W') then
	  nn = 0
	  d1dw = 0.0
          do l = js, je
	    do k = 1, nx
	      nn = nn + 1
	      d1dw(1:nz) = d1dw(1:nz) + w(k,l,1:nz)
	    enddo
	  enddo 
	  d1dw = d1dw / real(nn)    
!	        
        else if (trim(namev(i)) == 'PTIL' .or. trim(namev(i)) == 'QT' .or. trim(namev(i)) == 'QV'		 	&
	     .or. trim(namev(i)) == 'QC' .or. trim(namev(i)) == 'QR' .or. trim(namev(i)) == 'QI'			&
	     .or. trim(namev(i)) == 'NC' .or. trim(namev(i)) == 'NR' .or. trim(namev(i)) == 'NI'			&
	     .or. trim(namev(i)) == 'INDI' .or. trim(namev(i)) == 'INBCI' .or. trim(namev(i)) == 'INBAI'		&
	     .or. trim(namev(i)) == 'INDF' .or. trim(namev(i)) == 'INBCF' .or. trim(namev(i)) == 'INBAF' ) then
          print*, '    - Calculate and Write Resolved, Turbulent and precipitation fluxes: ', trim(namev(i))
	  namev2(mm+1) = namev(i)
	  nn = 0
	  d1d1 = 0.0
!
!  Net fluxes
!
     do l = js, je
	    do k = 1, nx
	      nn = nn + 1
	      d1d1(1:nz)   = d1d1(1:nz) + a(k,l,1:nz,i)
	      d1d2(1:nz,7*mm+1) = d1d2(1:nz,7*mm+1) + a(k,l,1:nz,i)*u(k,l,1:nz)
	      if (ny /= 1 .and. ivv == 1) d1d2(1:nz,7*mm+2) = d1d2(1:nz,7*mm+2) + a(k,l,1:nz,i)*v(k,l,1:nz)
	      d1d2(1:nz,7*mm+3) = d1d2(1:nz,7*mm+3) + a(k,l,1:nz,i)*w(k,l,1:nz)
	      ad1d2(k,l,1:nz,3*mm+1) = a(k,l,1:nz,i)*w(k,l,1:nz)
	    enddo
	  enddo
!
!  Turbulent fluxes
!
	  d1d1 = d1d1 / real(nn)
      do l = js, je
	    do k = 1, nx
	      d1d2(1:nz,7*mm+4) = d1d2(1:nz,7*mm+4) + (a(k,l,1:nz,i) - d1d1(1:nz))*(u(k,l,1:nz) - d1du(1:nz))
	      if (ny /= 1 .and. ivv == 1) d1d2(1:nz,7*mm+5) = d1d2(1:nz,7*mm+5) + (a(k,l,1:nz,i) - d1d1(1:nz))*(v(k,l,1:nz) - d1dv(1:nz))
	      d1d2(1:nz,7*mm+6) = d1d2(1:nz,7*mm+6) + (a(k,l,1:nz,i) - d1d1(1:nz))*(w(k,l,1:nz) - d1dw(1:nz))
	      ad1d2(k,l,1:nz,3*mm+2) = (a(k,l,1:nz,i) - d1d1(1:nz))*(w(k,l,1:nz) - d1dw(1:nz))
	    enddo
	  enddo
!
!  Precipitation fluxes
!
          if (trim(namev(i)) == 'QR' .or. trim(namev(i)) == 'QI' .or. trim(namev(i)) == 'NR' .or. trim(namev(i)) == 'NI'	&
	      .or. trim(namev(i)) == 'INDF' .or. trim(namev(i)) == 'INBCF' .or. trim(namev(i)) == 'INBAF') then
          do l = js, je
	    do k = 1, nx
	      vprecr = 0.; vpreci = 0.
	      if (iqr /= 0 .and. inr /= 0) vprecr(1:nz) = -572.212*(a(k,l,1:nz,inr)/a(k,l,1:nz,iqr))**(-0.3)
	      if (iqi /= 0 .and. ini /= 0) vpreci(1:nz) = -12.9743*(a(k,l,1:nz,ini)/a(k,l,1:nz,iqi))**(-0.2096)
	      if (trim(namev(i)) == 'QR' .or. trim(namev(i)) == 'NR') d1d2(1:nz,7*mm+7) = d1d2(1:nz,7*mm+7) + a(k,l,1:nz,i)*vprecr(1:nz)
	      if (trim(namev(i)) == 'QI' .or. trim(namev(i)) == 'NI') d1d2(1:nz,7*mm+7) = d1d2(1:nz,7*mm+7) + a(k,l,1:nz,i)*vpreci(1:nz)
	      if (trim(namev(i)) == 'INDF' .or. trim(namev(i)) == 'INBCF' .or. trim(namev(i)) == 'INBAF') 	&
	          d1d2(1:nz,7*mm+7) = d1d2(1:nz,7*mm+7) + a(k,l,1:nz,i)*vpreci(1:nz)
	      ad1d2(k,l,1:nz,3*mm+3) = a(k,l,1:nz,i)*vpreci(1:nz)
	    enddo
	  enddo
	  else
	    d1d2(1:nz,7*mm+7) = 0.
	    ad1d2(k,l,1:nz,3*mm+3) = 0.
	  endif
	  d1d2(:,7*mm+1:7*mm+7) = d1d2(:,7*mm+1:7*mm+7) / real(nn)
	  mm = mm + 1
	endif
      enddo
!   
!  Output
!  
      filename = 'fluxes.'//trim(time(j))//'.tec'
      open(13,FILE=trim(filename),FORM='formatted') 
      write(13,*) 'TITLE     = "fluxes "'//trim(time(j))
      write(13,*) 'VARIABLES = "x"'
      write(13,*) '"y"'
      write(13,*) '"z"'
      do i = 1, mm
        write(13,*) '"'//trim(namev2(i))//'"'
      enddo    
!      	  
      do k = 1, nz-1
        write(13,*) zz(k), real(d1d2(k,1:7*mm),4)
      enddo
      close (13)
!
      if (i3d == 1) then
      filename = 'vertflux3d.'//trim(time(j))//'.tec'
      open(14,FILE=trim(filename),FORM='formatted')   
      write(14,*) 'TITLE     = "fluxes 3D "'//trim(time(j))
      write(14,*) 'VARIABLES = "x"'
      write(14,*) '"y"'
      write(14,*) '"z"'
      do i = 1, mm
        if (trim(namev2(i)) == 'QI' .or. trim(namev2(i)) == 'NI'        					&
            .or. trim(namev2(i)) == 'INDF' .or. trim(namev2(i)) == 'INBCF' .or. trim(namev2(i)) == 'INBAF'	&
	    .or. trim(namev2(i)) == 'INDI' .or. trim(namev2(i)) == 'INBCI' .or. trim(namev2(i)) == 'INBAI') then
          write(14,*) '"'//trim(namev2(i))//'_NET"'
          write(14,*) '"'//trim(namev2(i))//'_TUR"'
          write(14,*) '"'//trim(namev2(i))//'_PRE"'
	endif
      enddo
      write(14,*) 'ZONE T="Zone"'
      write(14,*) 'I=',nx, ', J=', nny, ', K=', nz-1
      write(14,*) 'DATAPACKING=BLOCK'      
!
      write(14,101) xx(1:nx,js:je,1:nz-1)
      write(14,101) yy(1:nx,js:je,1:nz-1)
      write(14,101) zzz(1:nx,js:je,1:nz-1)      	  
      do i = 1, mm
        if (trim(namev2(i)) == 'QI' .or. trim(namev2(i)) == 'NI'      						&
            .or. trim(namev2(i)) == 'INDF' .or. trim(namev2(i)) == 'INBCF' .or. trim(namev2(i)) == 'INBAF'	&
	    .or. trim(namev2(i)) == 'INDI' .or. trim(namev2(i)) == 'INBCI' .or. trim(namev2(i)) == 'INBAI') then
          write(14,101) real(ad1d2(1:nx,js:je,1:nz-1,3*(i-1)+1),4)
	  write(14,101) real(ad1d2(1:nx,js:je,1:nz-1,3*(i-1)+2),4)
	  write(14,101) real(ad1d2(1:nx,js:je,1:nz-1,3*(i-1)+3),4)
	endif
      enddo
      close (14)
      endif
    endif
!
!  Read diagnostic files
!    
  if (j == 1 .and. idiagop == 1) then
    lfirst = 1
    idiag = 0
    sumid = 0
    do i = 1, nv
      ex5 = .FALSE.
      if (trim(namev(i)) /= 'Z0' .and. index(namev(i),'pro') == 0) then
        filename = dirname//trim(namev(i))//'.DIAG'
	if (trim(namev(i)) == 'V' .and. ny == 1) filename = dirname//'NOV.DIAG'
!
      inquire (FILE=trim(filename), EXIST=ex5)
      if (ex5) then
        idiag(i) = 1
      endif
      sumid = sumid + idiag(i)
!
      if (ex5) then
      ifile = 13+i
      open(ifile,FILE=trim(filename),FORM=form,STATUS='OLD')
      if (lfirst == 1) then
        ierr = 0
        itt = 0
        do while (ierr == 0) 
          read(ifile,'(1x,a1)',iostat=ierr) car1
          if (car1 == '#') then
	    read(ifile,*) tt
	    if (tt >= stime) itt = itt + 1
	  endif
        enddo 
	lfirst = 0     
      endif
      endif
      endif
    enddo
!
    allocate (a2dd(1:itt,1:nz,1:15,1:sumid+3))
!
    k = 0
    print*, ''
    do i = 1, nv
      if (idiag(i) == 1) then
      k = k + 1
      ifile = 13+i
      filename = dirname//trim(namev(i))//'.DIAG'
      print*, 'Open and read Diagnostic file ', trim(filename)
      rewind (ifile)
!
      if (trim(namev(i)) == 'QI') iqid = k
!
      do m = 1, itt
        do n = 1, 15
          a2dd(m,:,n,sumid+2) = zzz(2,js,:)
	enddo
        a2dd(m,:,:,sumid+3) = 0.
	read(ifile,'(1x,a1)') car1
98      read(ifile,'(e15.7)') tt
	if (tt >= stime) then
	  a2dd(m,:,:,sumid+1) = tt
          if (iform == 1) then
	    do n = 1, 15
	      read(ifile,*)
    	      read(ifile,100) a2dd(m,:,n,k)   
	    enddo   
          else
	    do n = 1, 15
	      read(ifile,*)
  	      read(ifile) a2dd(m,:,n,k)
	    enddo    
          endif
	else
	  read(ifile,'(1x,a1)') car1
	  do while (car1 /= '#')
	  read(ifile,'(1x,a1)') car1
	  enddo
	  goto 98
	endif
      enddo
      close (ifile)      
      endif
    enddo 
!
!  Write Diagnostic files
!
    if (sumid /= 0) then
      filename = 'diagnostics.tec'
      open(13,FILE=trim(filename),FORM='formatted')  
!      
      write(13,*) 'TITLE     = "diagnostics"'
      write(13,*) 'VARIABLES = "x"'
      write(13,*) '"y"'
      write(13,*) '"z"'
      do i = 1, nv
        if (idiag(i) == 1) then
	  do k = 1, 15
   	    write(car2,'(i2)') k
            if (trim(namev(i)) /= 'Z0' .and. trim(namev(i)) /= 'K' .and. index(namev(i),'pro') == 0)    &
              write(13,*) '"'//trim(namev(i))//'.'//trim(car2)//'"'
	  enddo
	endif
      enddo
      write(13,*) 'ZONE T="Zone"'
      write(13,*) 'I=', itt, ', J=', 1, ', K=', nz-1
      write(13,*) 'DATAPACKING=BLOCK'
!
      write(13,101) real(a2dd(1:itt,1:nz-1,1,sumid+1),4)
      write(13,101) real(a2dd(1:itt,1:nz-1,1,sumid+3),4)
      write(13,101) real(a2dd(1:itt,1:nz-1,1,sumid+2),4)
      do i = 1, sumid           
        do k = 1, 15
          write(13,101) real(a2dd(1:itt,1:nz-1,k,i),4)
	enddo
      enddo
      close (13)
    endif
!
  endif
!
!  Read 2d profiles (z-t)
!
  if (j == 1) then
!
  print*, ''
  ex4 = .FALSE.
  iprof = 0
  lfirst = 1
  ispro = 2
  do i = 1, nv
    if (trim(namev(i)) /= 'Z0') then
      filename = dirname//trim(namev(i))//'.0000'
      if (trim(namev(i)) == 'T.pro') item = i
      if (trim(namev(i)) == 'QC.pro') iqc = i
      if (trim(namev(i)) == 'NC.pro') inc = i
      if (trim(namev(i)) == 'QR.pro') iqr = i
      if (trim(namev(i)) == 'NR.pro') inr = i
      if (trim(namev(i)) == 'QI.pro') iqi = i
      if (trim(namev(i)) == 'NI.pro') ini = i
      if (trim(namev(i)) == 'INDF.pro') indf = i
      if (trim(namev(i)) == 'INBCF.pro') inbcf = i
      if (trim(namev(i)) == 'RADDT.pro') irad = i
    else
      cycle
    endif
    inquire (FILE=trim(filename), EXIST=ex4)    
!  
    if (index(namev(i),'pro') == 0) ispro = ispro + 1
!
    if (ex4) then
      ifile = 12+i
      print*, 'Open and read in 2d profile file ', trim(filename)
      open(ifile,FILE=trim(filename),FORM=form,STATUS='OLD')
      if (lfirst == 1) then
        ierr = 0
	itt = 0
        do while (ierr == 0) 
	  read(ifile,'(1x,a1)',iostat=ierr) car1
	  if (car1 == '#') then
  	    read(ifile,*) tt
	    if (tt >= stime) itt = itt + 1
	  endif
        enddo
!
	rewind (ifile)
        allocate(a2d(1:itt,1:nz,1:nv+3))
        lfirst = 0
      endif
!
      iprof = iprof+1
      do m = 1, itt
        a2d(m,:,nv+2) = zzz(2,js,:)
        a2d(m,:,nv+3) = 0.
	read(ifile,'(1x,a1)') car1
99      read(ifile,'(e15.7)') tt
	if (tt >= stime) then
	  a2d(m,:,nv+1) = tt
          if (iform == 1) then
  	    read(ifile,100) a2d(m,:,i)      
          else
  	    read(ifile) a2d(m,:,i)
          endif
	else
	  read(ifile,'(1x,a1)') car1
	  do while (car1 /= '#')
	  read(ifile,'(1x,a1)') car1
	  enddo
	  goto 99	
	endif
      enddo
      print*, '	 - Vertical profile in 2D (z - t)'    
      close (ifile)
    endif
  enddo
!
!  Calculate derived time series
!
  if (iseries == 1 .and. iprof2d == 1) then
    allocate (series(1:itt,1:10))
    series = 0.
    series(:,8) = 1000.
    do m = 1, itt
      sum = 0.
      do k = 2, nz-1
        if (a2d(m,k,inc) > 10.) then
	  series(m,1) = series(m,1) + a2d(m,k,inc)
	  series(m,2) = series(m,2) + a2d(m,k,iqc)
	  series(m,3) = series(m,3) + a2d(m,k,ini)
	  series(m,4) = series(m,4) + a2d(m,k,iqi)
	  series(m,5) = series(m,5) + a2d(m,k,indf)
	  series(m,6) = series(m,6) + a2d(m,k,inbcf)
	  sum = sum + 1.
	endif
	if (a2d(m,k,irad) < series(m,7)) series(m,7) = a2d(m,k,irad)
	if (a2d(m,k,item) < series(m,8) .and. a2d(m,k,nv+2) < 1200.) series(m,8) = a2d(m,k,item)
	if (a2dd(m,k,12,iqid) > series(m,9) .and. a2d(m,k,nv+2) < 1200.) series(m,9) = a2dd(m,k,12,iqid)
	if ((a2dd(m,k,13,iqid) > series(m,10) .and. a2dd(m,k,13,iqid) > a2dd(m,k+1,13,iqid)) .or.		&
	    (a2dd(m,k,13,iqid) < series(m,10) .and. a2dd(m,k,13,iqid) < a2dd(m,k+1,13,iqid)) .and. 		&
	    a2d(m,k,nv+2) < 1200.) series(m,10) = a2dd(m,k,13,iqid)
      enddo
      series(m,1:6) = series(m,1:6)/max(sum,1.)
    enddo
  endif
!
!  Write vertical profiles in 2D (z - t)
!
    if (iprof2d == 1 .and. ex4) then
      filename = 'profiles.2d.tec'
      open(12,FILE=trim(filename),FORM='formatted')  
!      
      write(12,*) 'TITLE     = "profiles 2d"'
      write(12,*) 'VARIABLES = "x"'
      write(12,*) '"y"'
      write(12,*) '"z"'
      do i = 1, nv
        if (trim(namev(i)) /= 'Z0' .and. trim(namev(i)) /= 'K' .and. index(namev(i),'pro') /= 0) then
           write(12,*) '"'//trim(namev(i))//'"'
        endif
      if (trim(namev(i)) == 'K')     &
        write(12,*) '"Ksgs"'
      enddo
      write(12,*) 'ZONE T="Zone"'
      write(12,*) 'I=', itt, ', J=', 1, ', K=', nz-1
      write(12,*) 'DATAPACKING=BLOCK'
!
      write(12,101) real(a2d(1:itt,1:nz-1,nv+1),4)
      write(12,101) real(a2d(1:itt,1:nz-1,nv+3),4)
      write(12,101) real(a2d(1:itt,1:nz-1,nv+2),4)
      do i = ispro, nv              
          write(12,101) real(a2d(1:itt,1:nz-1,i),4)
      enddo 
!
!  Write derived time series
!
      if (iseries == 1) then
        filename = 'time_series_derived.dat'
        open(11,FILE=trim(filename),FORM='formatted')  
        do m = 1, itt
          write(11,*) a2d(m,1,nv+1), series(m,:)
        enddo 
        close(11)  
      endif
!
!  Write regular profiles
!    
    else if (iprof2d == 0 .and. ex4) then
      allocate (time2(1:itt))
      do m = 1, itt
        write(time2(m),'(i5)') floor(a2d(m,1,nv+1))
        filename = 'profiles.'//trim(adjustl(time2(m)))//'.tec'
        open(11,FILE=trim(filename),FORM='formatted')  
        do i = 1, nz-1
          write(11,102) real(zz(i),4), real(a2d(m,i,ispro:nv))
        enddo 
        close(11)  
      enddo  
    endif
!
  endif
!
!  Write vertical profiles
!
  if (.not.ex4) then
    filename = 'profiles.'//trim(time(j))//'.tec'
    open(11,FILE=trim(filename),FORM='formatted')  
    do i = 1, nz-1
      write(11,102) real(zz(i),4), real(a1d(i,1:iprof))
    enddo 
    close(11)
  endif
!
100  format(6e15.7)
101  format(3(e15.7,3x))  
102  format(e15.7,5x, 38(e15.7,3x))  
103  format(38(e15.7,3x))
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
  inquire(UNIT=ifile, OPENED=ll)
  if (ll) close(ifile)
  inquire(UNIT=10, OPENED=ll)
  if (ll) close(10)
  inquire(UNIT=12, OPENED=ll)
  if (ll) close(12)
!
  enddo
!
end program
