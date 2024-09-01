
#include "ctrparam.h"
  
! ================================================================
!
!  Purpose:
!      A simple online 2D cloud tracking algorithm
!
!  Author
!      Julien Savre
!      MIM, Ludwig-Maximilian Universität, Munich
!
! ================================================================

!      =================================================                            
  subroutine cloud_tracking ( dens, wind, thermo, hydromtr, lout )
!      =================================================                            
  
  USE gridno
  USE shared_data
  USE shared_wind
  USE shared_hydro
  USE shared_thermo
  USE shared_pressure
  USE netcdfmod
  USE averages

  USE searchmod
  USE shared_clouds

#ifdef SPMD
  USE mpicomm
#endif

  IMPLICIT NONE

  integer :: unit1, unit2, unit3, io
  integer :: n, i, j, k
  logical :: ex, lout, lfirst=.true.
  character(len = 100) :: outputname, filename
  character(len = 4) :: car4
  character(len = 1) :: car1
!
  real :: zz, tav, nav, navb, nupav, nupavb, nvar, nvarb, nupvar, nupvarb
  real :: mav, mavb, mvar, mvarb, mupav, mupavb, mupvar, mupvarb, msummax, mupsummax
  real :: dav, dvar, dupav, dupvar
  real, dimension(npdf) :: pdf, pdfm, pdfupm, pdfd, pdfupd, pdft
  real, dimension(npdf/2,npdf/2) :: pdfda, pdfdd
  real, dimension(npdf,2) :: pdfdq
!
  real, dimension (ip_start:ip_end,jp_start:jp_end,1:nz) :: dens
  type (winds) :: wind
  type (thermodyn) :: thermo
  type (hydrometeor), dimension (1:nhydro)	:: hydromtr   
!
  real, dimension (2) :: tt
  real, dimension (1:nz)  :: bav
  real, dimension (ip_start:ip_end,jp_start:jp_end) :: mom, buo, liq, lwp, ctop
  real, dimension (1:nx,1:ny) :: mom_all, buo_all, liq_all, lwp_all, ctop_all, cl_id, up_id
  integer, dimension (ip_start:ip_end,jp_start:jp_end) :: if1, if2
  integer, dimension (1:nx,1:ny) :: if1_all, if2_all
!
  type (cloud_type), dimension(ncl) :: clouds_old

  namelist /clt/object,nt,dxl,dyl,h_clt,qcmin,wmin,bmin,lwpmin,ctmin,outputname,nc_min,nu_min,nblock,time_clt,split_clt,stat_clt,with_up
!
  save lfirst, tt
!
!----------------------------------------------------------!
!                 First, read input file                   !
!----------------------------------------------------------!
!
  if (lfirst) then
    allocate( clouds(1:ncl), updrafts(1:ncl,1:nup) )
!
    filename = 'namelist'
    open(3,file=filename(1:len_trim(filename)),status='old')
    read(3,nml=clt)
    close(3)
  endif
!
  allocate( nupdrafts(1:ncl) )
!
  NXO = nx-5
  is = 4
  ie = nx-2
#ifdef MODEL_3D
  NYO = ny-5
  js = 4
  je = ny-2
#else
  NYO = 1
  js = 1
  je = 1
#endif
  NZO = nz
!
  tt(2) = time
  if (lfirst) tt(1) = time
!
!----------------------------------------------------------!
!              Initialize clouds data type                 !
!----------------------------------------------------------!
!
    if (lfirst) then
      do i = 1, ncl
        clouds_old(i) = 0.
      enddo
      clouds%age=0.
    else
      do i = 1, ncl
        clouds_old(i) = clouds(i)
      enddo
    endif
!
    do i = 1, ncl
      clouds(i) = 0.
      do j = 1, nup
        updrafts(i,j) = 0.
      enddo
    enddo
!
!----------------------------------------------------------!
!                    Initialize slice                      !
!----------------------------------------------------------!
!
!  Projecting clouds on a 2D surface
!
    if1=0
    if2=0
    cl_id = 0
    up_id = 0
    lwp = 0.
    ctop = 0.
!
    call horav (thermo%buoy, bav)
!
    zz = 0.5*dz*fdz0(1)
    do k=1,nz-1
      if ( zz >= h_clt .and. zz-dz*fdz0(k) < h_clt ) then
        mom = 0.5*dens(:,:,k)*(wind(:,:,k+1)%w + wind(:,:,k)%w)
        liq = hydromtr(:,:,k,drop)%q + hydromtr(:,:,k,rain)%q
        buo = thermo%buoy(:,:,k) - bav(k)
      endif
      where ( hydromtr(:,:,k,drop)%q + hydromtr(:,:,k,rain)%q > qcmin ) 
        lwp = lwp + dens(:,:,k)*(hydromtr(:,:,k,drop)%q + hydromtr(:,:,k,rain)%q)*dz*fdz0(k)
        ctop = zz
      end where
!
      if ( (trim(object) == 'cloud_base' .or. trim(object) == 'cold_pool') ) then  		  
!
        if ( zz >= h_clt ) then
	  if (trim(object) == 'cloud_base') then
     	    where ( liq > qcmin ) if1 = 1
   	    where ( mom >= wmin .and. buo >= bmin ) if2 = 1
   	  else if (trim(object) == 'cold_pool') then
     	    where ( buo < bmin ) if1 = 1
	  endif
	  exit
	endif
!
      else if ( trim(object) == 'cloud_proj' ) then
!
	if (k == nz-1) then
	  where (lwp > lwpmin .and. ctop >= ctmin) if1 = 1
	  where (lwp > lwpmin .and. ctop < ctmin) if2 = 1
	endif
!
      else 
!
        print*, 'Error, object type not recognized'
	stop
!
      endif
!
      if (k < nz) zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))
    enddo
!
!  Gather on master proc
!
#ifdef SPMD
  call collect ( if1, if1_all )
  call collect ( if2, if2_all )
  call collect ( mom, mom_all )
  call collect ( buo, buo_all )
  call collect ( liq, liq_all )
  call collect ( lwp, lwp_all )
  call collect ( ctop, ctop_all )
#else
  if1_all = if1
  if2_all = if2
  mom_all = mom
  buo_all = buo
  liq_all = liq
  lwp_all = lwp
  ctop_all = ctop
#endif
!
!----------------------------------------------------------!
!	         Clouds and updrafts search                !
!----------------------------------------------------------!
!
  if ( mypid == 0 ) then
    call object_search ( if1_all, cl_id, clouds, .true. )
!
    if (with_up) call subobject_search ( if1_all, if2_all, up_id, clouds, updrafts )
!
    if (with_up .and. split_clt) call object_split ( cl_id, clouds, updrafts )
!
    if (.not.lfirst .and. time_clt) call time_tracking ( 2, tt, clouds, clouds_old, cl_id )
  endif
!
!----------------------------------------------------------!
!		   Diagnostics and outputs                 !
!----------------------------------------------------------!
!
!  Output slice
!
  if ( mypid == 0 .and. lout ) call write_tmp_slice ( 0, h_clt, mom_all, buo_all, liq_all, lwp_all, ctop_all, cl_id, up_id )
!
  if ( mypid == 0 .and. lout .and. stat_clt ) then
!
!  Open output files
!
    write(car4,'(i4)') floor(h_clt)
!
    unit1=203
    unit2=204
    unit2=205
!
    filename = './CL_TS' // adjustl(adjustr(car4))
    inquire (FILE=filename, OPENED=ex)
    if (.not.ex) then
      open(unit1,file=filename(1:len_trim(filename)), 	   &
  	     form='formatted',				   &
  	     access='SEQUENTIAL', 			   &
  	     position='APPEND',			   &
  	     status='unknown')
    endif
!
    filename = './CL_PDF' // adjustl(adjustr(car4))
    inquire (FILE=filename, OPENED=ex)
    if (.not.ex) then
      open(unit2,file=filename(1:len_trim(filename)), 	   &
  	     form='formatted',				   &
  	     access='SEQUENTIAL', 			   &
  	     position='REWIND',				   &
  	     status='unknown')
    endif
!
    filename = './CL_JPDF' // adjustl(adjustr(car4))
    inquire (FILE=filename, OPENED=ex)
    if (.not.ex) then
      open(unit3,file=filename(1:len_trim(filename)), 	   &
  	     form='formatted',				   &
  	     access='SEQUENTIAL', 			   &
  	     position='REWIND',				   &
  	     status='unknown')
    endif
!
!  Calculate diagnostics
!
    call cloud_diagnostics ( mom_all, clouds, updrafts )
!
    call cloud_stats ( clouds, clouds_old, updrafts, tav, nav, navb, 			&
    		       nupav, nupavb, nvar, nvarb, nupvar, nupvarb, pdf, 		&
		       mav, mavb, mvar, mvarb, pdfm, mupav, mupavb, mupvar, mupvarb, 	&
		       msummax, mupsummax, pdfupm, dav, dvar, 				&
		       pdfd, dupav, dupvar, pdfupd, pdft, pdfdq, pdfda, pdfdd )
!
!  Outputs
!
    call cloud_output ( unit1, unit2, unit3, time, tav, maxval(clouds%size),		&
    			nav, nvar, navb, nvarb, nupav, nupvar, 				&
			nupavb, nupvarb, mav, mvar, mavb, mvarb, 			&
			mupav, mupvar, mupavb, mupvarb, msummax, mupsummax, 		&
			dav, dvar, dupav, dupvar, pdf, pdfm, pdfupm, pdfd, pdfupd, 	&
			pdft, pdfdq, pdfda, pdfdd )
!
    close (unit1)
    close (unit2)
    close (unit3)
!
  endif
!
  if (lfirst) lfirst=.false.
!
  deallocate( nupdrafts )
!
  tt(1) = tt(2)
!
!----------------------------------------------------------!
!  
  end subroutine cloud_tracking
