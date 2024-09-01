  program particule_reader
!  
!----------------------------------------------------------!
!  
    IMPLICIT NONE
!    
    integer :: npart, ntime
    integer :: np, nt, k, l, p1, eof
    character(len=100) :: car100
!    
    integer :: i1, i2, mask
    integer, dimension(:,:), allocatable :: mask_t
    real :: time, x, y, z, u, v, w, pt, qt, qc, rh, buo
    real, dimension(:,:), allocatable :: time_t, x_t, y_t, z_t, u_t, v_t, w_t, pt_t, qt_t, qc_t, rh_t, buo_t
    real, dimension(:,:), allocatable :: h_ent, t_ent, h_det, t_det
!
    integer :: ii, ipdf(50,3)
    real :: dpdf
    real, dimension(50,3) :: pdfp, pdfq
!
    real, parameter :: qcmin=1.e-6, wmin=0.1, zmin=3000.
!  
!----------------------------------------------------------!
!   	       Opening and reading particles		   !
!----------------------------------------------------------!
!
  open(UNIT=100,FILE='./OUTPUT/PARCELS',FORM='formatted',STATUS='old')
!
  eof = 0
  ntime = 0
  do 
    ntime = ntime + 1
    read(100,*,iostat=eof) npart
    if ( eof /= 0 ) exit
    do np = 1, npart
      read(100,101) time, i1, i2, mask, x, y, z, u, v, w, pt, qt, qc, rh, buo
    enddo
  enddo
  ntime = ntime-1
  rewind(100)
!
  allocate( x_t(1:npart,1:ntime), y_t(1:npart,1:ntime), z_t(1:npart,1:ntime) ) 
  allocate( u_t(1:npart,1:ntime), v_t(1:npart,1:ntime), w_t(1:npart,1:ntime) )
  allocate( qt_t(1:npart,1:ntime), qc_t(1:npart,1:ntime), rh_t(1:npart,1:ntime) )
  allocate( pt_t(1:npart,1:ntime), buo_t(1:npart,1:ntime) )
  allocate( time_t(1:npart,1:ntime), mask_t(1:npart,1:ntime) )
! 
  do nt = 1, ntime
    read(100,*) npart
    do np = 1, npart
      read(100,101) time_t(np,nt), i1, i2, mask_t(np,nt), x_t(np,nt), y_t(np,nt), z_t(np,nt), u_t(np,nt), v_t(np,nt), w_t(np,nt),  &
      	pt_t(np,nt), qt_t(np,nt), qc_t(np,nt), rh_t(np,nt), buo_t(np,nt)
    enddo
  enddo
! 
  close(100)
!
  print*,  npart, 'particles over', ntime, 'time steps read'
!  
!----------------------------------------------------------!
!   	              Post processing			   !
!----------------------------------------------------------!
!
  allocate( h_ent(1:npart,1:50), t_ent(1:npart,1:50), h_det(1:npart,1:50), t_det(1:npart,1:50) )
!
  h_ent = 0.; t_ent = 0.; h_det = 0.; t_det = 0.
!
  ipdf = 0
  pdfp = 0.
  pdfq = 0.
  dpdf = (12000. - zmin)/50.
  do np = 1, npart
    k = 0
    l = 0
    do nt = 2, ntime
!
      if ( z_t(np,nt) > zmin ) then
        if ( qc_t(np,nt-1) < qcmin .and. qc_t(np,nt) >= qcmin ) then 
	  k = k + 1
	  h_ent(np,k) = z_t(np,nt)
	  t_ent(np,k) = time_t(np,nt)
        else if ( qc_t(np,nt-1) >= qcmin .and. qc_t(np,nt) < qcmin ) then 
	  l = l + 1
	  h_det(np,l) = z_t(np,nt)
	  t_det(np,l) = time_t(np,nt)	       
        endif
      endif
!
      if ( z_t(np,nt) > zmin ) then
        if ( qc_t(np,nt-1) < qcmin .and. qc_t(np,nt) >= qcmin ) then 
	  ii = floor((z_t(np,nt) - zmin)/dpdf) + 1
!
	  if (time_t(np,nt) >= time_t(np,ntime)/4. .and. time_t(np,nt) < time_t(np,ntime)*2./4.) then
	    pdfp(ii,1) = pdfp(ii,1) + pt_t(np,nt-1)*(1. + 0.601*(qt_t(np,nt-1) - qc_t(np,nt-1)))
	    pdfq(ii,1) = pdfq(ii,1) + rh_t(np,nt-1)
	    ipdf(ii,1) = ipdf(ii,1) + 1
	  else if (time_t(np,nt) >= time_t(np,ntime)*2./4. .and. time_t(np,nt) < time_t(np,ntime)*3./4.) then
	    pdfp(ii,2) = pdfp(ii,2) + pt_t(np,nt-1)*(1. + 0.601*(qt_t(np,nt-1) - qc_t(np,nt-1)))
	    pdfq(ii,2) = pdfq(ii,2) + rh_t(np,nt-1)
	    ipdf(ii,2) = ipdf(ii,2) + 1
	  else if (time_t(np,nt) >= time_t(np,ntime)*3./4.) then
	    pdfp(ii,3) = pdfp(ii,3) + pt_t(np,nt-1)*(1. + 0.601*(qt_t(np,nt-1) - qc_t(np,nt-1)))
	    pdfq(ii,3) = pdfq(ii,3) + rh_t(np,nt-1)
	    ipdf(ii,3) = ipdf(ii,3) + 1
	  endif
        endif
      endif
!      
    enddo
  enddo
!
!  Print results
!
  do
!	
    print*, 
    print*, 'Choose a particle number, or exit (q)'
!
    read*, car100
    if (trim(car100) == 'q') exit
    read(car100,'(i7)')p1
    if (p1 < 0 .or. p1 > npart) cycle
!
    open(UNIT=99,FILE='./OUTPUT/PARCEL_SINGLE',FORM='formatted',STATUS='unknown')
    do nt = 1, ntime
      write(99,101) time_t(p1,nt), i1, i2, mask_t(p1,nt), x_t(p1,nt), y_t(p1,nt), z_t(p1,nt), u_t(p1,nt), v_t(p1,nt), w_t(p1,nt),  &
      	pt_t(p1,nt), qt_t(p1,nt), qc_t(p1,nt), rh_t(p1,nt), buo_t(p1,nt)
    enddo
    close(99)
!
    do k = 1, 50
      write(15,*) zmin+real(k-1)*dpdf, pdfp(k,1)/real(ipdf(k,1)), pdfp(k,2)/real(ipdf(k,2)), pdfp(k,3)/real(ipdf(k,3)), &
      	pdfq(k,1)/real(ipdf(k,1)), pdfq(k,2)/real(ipdf(k,2)), pdfq(k,3)/real(ipdf(k,3))
    enddo
!
  enddo
!
  deallocate( h_ent, t_ent, h_det, t_det )
!  
!----------------------------------------------------------!
!   	      		 Closing			   !
!----------------------------------------------------------!
!
  deallocate( time_t,x_t,y_t,z_t,u_t,v_t,w_t )
  deallocate( mask_t,qt_t,qc_t,rh_t,pt_t,buo_t )
!
101  format (2x,f16.8,3i16,11f16.8)
!
!----------------------------------------------------------!
!  
  end program
