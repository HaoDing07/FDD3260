
!      =================================================                            
  program cloud_tracking
!      =================================================                            

  USE sharedmod
  USE cloudsmod
  USE netcdfmod
  USE searchmod

  IMPLICIT NONE

  integer :: unit1, unit2, unit3, io
  integer :: n, i, j, k, t, nenv, nclo
  logical :: ex, lfirst=.true.
  character(len = 100) :: outputname, filename
  character(len = 4) :: car4
  character(len = 1) :: car1
!
  real :: nt, dt, rho0
  real :: x0, tav, nav, navb, nupav, nupavb, nvar, nvarb, nupvar, nupvarb, qe, qc
  real :: mav, mavb, mvar, mvarb, mupav, mupavb, mupvar, mupvarb, msummax, mupsummax
  real :: dav, dvar, dupav, dupvar
  real, dimension(npdf) :: pdf, pdfm, pdfupm, pdfd, pdfupd, pdft
  real, dimension(npdf/2,npdf/2) :: pdfda, pdfdd
  real, dimension(npdf,2) :: pdfdq
!
  type (cloud_type), dimension(ncl) :: clouds, clouds_old
  type (cloud_type), dimension(ncl,nup) :: updrafts
!
  real, dimension(:,:), allocatable :: mom_all, qt_all, buo_all, bef_all, mse_all, liq_all, lwp_all, ctop_all, ent_all, det_all
  real, dimension(:,:,:), allocatable :: dq_all
  real, dimension(:,:), allocatable :: cl_id, cl_rad, cl_mfl, cl_mse, cl_buo, cl_ent, cl_det, up_id
  integer, dimension(:,:), allocatable :: if1, if2

  namelist /clt/object,nt,dxl,dyl,h_clt,rho0,qcmin,wmin,bmin,lwpmin,ctmin,outputname,nc_min,nu_min,nblock,time_clt,stat_clt,split_clt,with_up,with_halo
!
  save lfirst
!
!----------------------------------------------------------!
!                 First, read input file                   !
!----------------------------------------------------------!
!
!  Option namelist
!
  filename = 'namelist'
  open(3,file=filename(1:len_trim(filename)),status='old')
  read(3,nml=clt)
  close(3)
!
  call read_dim (outputname)
!
!  Local allocations
!
  allocate( mom_all(1:NXO,1:NYO), buo_all(1:NXO,1:NYO), bef_all(1:NXO,1:NYO),  &
  	mse_all(1:NXO,1:NYO), liq_all(1:NXO,1:NYO), lwp_all(1:NXO,1:NYO),      &
	ctop_all(1:NXO,1:NYO), ent_all(1:NXO,1:NYO), det_all(1:NXO,1:NYO),     &
	dq_all(1:NXO,1:NYO,5), qt_all(1:NXO,1:NYO), if1(1:NXO,1:NYO), if2(1:NXO,1:NYO) )
  allocate( cl_id(1:NXO,1:NYO), cl_rad(1:NXO,1:NYO), cl_mfl(1:NXO,1:NYO),      &
        cl_mse(1:NXO,1:NYO), cl_buo(1:NXO,1:NYO), cl_ent(1:NXO,1:NYO), 	       &
	cl_det(1:NXO,1:NYO), up_id(1:NXO,1:NYO) )
  allocate( nupdrafts(1:ncl) )
!
  is = 1
  ie = NXO
  js = 1
  je = NYO
!
  dt = 60.
!
!----------------------------------------------------------!
!              Initialize clouds data type                 !
!----------------------------------------------------------!
!
  do t = 1, NTO
!
    if ( nt <= 0. .or. (time(t) > nt-dt .and. time(t) < nt+dt) ) then
!
    print*, "Extracting cloud objects at time",  time(t)
!
!  Clouds
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
!
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
    qe = 0.
    qc = 0.
    cl_id = 0.
    cl_rad= 0.
    up_id = 0.
    nenv = 0
    nclo = 0
!
    do k=1,NZO
      x0 = sum(state(1:NXO,1:NYO,k,t)%buoy) / real(NXO*NYO)
      do j=1,NYO
        do i=1,NXO
!
	  if (rho0 /= 0.) then
	    mom_all(i,j) = rho0*wind(i,j,k,t)%w
	  else
            mom_all(i,j) = 0.5*state(i,j,k,t)%dens*(wind(i,j,k+1,t)%w + wind(i,j,k,t)%w)
	  endif
          liq_all(i,j) = state(i,j,k,t)%qc + state(i,j,k,t)%qr
          qt_all(i,j)  = state(i,j,k,t)%qv + state(i,j,k,t)%qc + state(i,j,k,t)%qr
          buo_all(i,j) = state(i,j,k,t)%buoy - x0
          bef_all(i,j) = state(i,j,k,t)%beff
          mse_all(i,j)  = state(i,j,k,t)%mse
          ent_all(i,j)  = state(i,j,k,t)%ent
          det_all(i,j)  = state(i,j,k,t)%det
          lwp_all(i,j)  = state(i,j,k,t)%lwp
          ctop_all(i,j) = state(i,j,k,t)%ctop
          dq_all(i,j,1) = tend(i,j,k,t)%dqcon
          dq_all(i,j,2) = tend(i,j,k,t)%dqadv
          dq_all(i,j,3) = tend(i,j,k,t)%dqdif
          dq_all(i,j,4) = tend(i,j,k,t)%dqfor
          dq_all(i,j,5) = tend(i,j,k,t)%dqdt
!
          if ( (trim(object) == 'cloud_base' .or. trim(object) == 'cold_pool') ) then  		  
!
            if ( zz(k) >= h_clt .or. NZO == 1 ) then
	      if (trim(object) == 'cloud_base') then
!     	        if ( liq_all(i,j) > qcmin ) if1(i,j) = 1
     	        if ( mom_all(i,j) >= wmin .and. liq_all(i,j) > qcmin ) if1(i,j) = 1
   	        if ( mom_all(i,j) >= wmin .and. buo_all(i,j) >= bmin ) if2(i,j) = 1
		if ( if1(i,j) < 0.5 ) qe = qe + qt_all(i,j)
		if ( if1(i,j) < 0.5 ) nenv = nenv + 1
		if ( if1(i,j) > 0.5 ) qc = qc + qt_all(i,j)
		if ( if1(i,j) > 0.5 ) nclo = nclo + 1
   	      else if (trim(object) == 'cold_pool') then
     	        if ( buo_all(i,j) < bmin ) if1(i,j) = 1
	      endif
	    endif
!
          else if ( trim(object) == 'cloud_proj' ) then
!
	    if (k == NZO) then
	      if (lwp_all(i,j) > lwpmin .and. ctop_all(i,j) >= ctmin) if1(i,j) = 1
	      if (lwp_all(i,j) > lwpmin .and. ctop_all(i,j) < ctmin) if2(i,j) = 1
	    endif
!
          else 
!
            print*, 'Error, object type not recognized'
	    stop
!
          endif
!
        enddo
      enddo
    enddo
    qe = qe / real(nenv)
    qc = qc / real(nclo)
!
!----------------------------------------------------------!
!	         Clouds and updrafts search                !
!----------------------------------------------------------!
!
    call object_search ( if1, cl_id, cl_rad, clouds )
!
    if (with_halo) call find_halos ( cl_id, clouds )
!
    if (with_up) call subobject_search ( if1, if2, up_id, clouds, updrafts )
!
    if (with_up .and. split_clt) call object_split ( cl_id, clouds, updrafts )
!
    if (.not.lfirst .and. time_clt) call time_tracking ( t, time, clouds, clouds_old, cl_id )
!
!----------------------------------------------------------!
!		        Diagnostics                        !
!----------------------------------------------------------!
!
!  Calculate diagnostics
!
    call cloud_diagnostics ( mom_all, qt_all, mse_all, buo_all, bef_all, ent_all, det_all, lwp_all, ctop_all, dq_all, qe, clouds, updrafts )
!
    do n = 1, ncl
      where ( cl_id == clouds(n)%id ) 
        cl_mfl = clouds(n)%mfl
	cl_mse = clouds(n)%mse
	cl_buo = clouds(n)%buo
	cl_ent = clouds(n)%ent
	cl_det = clouds(n)%det
      end where
    enddo
!
!----------------------------------------------------------!
!		          Outputs                          !
!----------------------------------------------------------!
!
!  Output slice
!
    call write_slice ( t, h_clt, mom_all, buo_all, mse_all, liq_all, lwp_all, 	&
    		       ctop_all, cl_id, cl_rad, cl_mfl, cl_mse, cl_buo, cl_ent, cl_det, up_id )
!
!  Output cloud scatter plots
!
    open(202, file='cloud_scatter.dat', status='unknown')
    do n = 1, ncl
      if ( clouds(n)%size > nc_min )	&
        write(202,*) sqrt(dxl*dyl*real(clouds(n)%size)), clouds(n)%sf, clouds(n)%mfl, clouds(n)%mom, clouds(n)%mse, clouds(n)%buo, clouds(n)%bef, clouds(n)%ent/clouds(n)%mfl, -clouds(n)%det/clouds(n)%mfl, clouds(n)%lwp, clouds(n)%ctop, max(clouds(n)%dmdz,0.)/(qe*clouds(n)%mfl), -min(clouds(n)%dmdz,0.)/(qc*clouds(n)%mfl)
    enddo
    close(202)
!
!  Calculate and output statistics
!
    if (stat_clt) then
!
      unit1=203
      unit2=204
      unit2=205
!
      write(car4,'(i4)') floor(h_clt)
!
      filename = './CL_TS' // adjustl(adjustr(car4))
      inquire (FILE=filename, OPENED=ex)
      if (.not.ex) then
        open(unit1,file=filename(1:len_trim(filename)), 	&
  	       form='formatted',				&
  	       access='SEQUENTIAL', 			        &
  	       position='APPEND',			        &
  	       status='unknown')
      endif
!
      filename = './CL_PDF' // adjustl(adjustr(car4))
      inquire (FILE=filename, OPENED=ex)
      if (.not.ex) then
        open(unit2,file=filename(1:len_trim(filename)), 	&
  	       form='formatted',				&
  	       access='SEQUENTIAL', 			   	&
  	       position='REWIND',				&
  	       status='unknown')
      endif
!
      filename = './CL_JPDF' // adjustl(adjustr(car4))
      inquire (FILE=filename, OPENED=ex)
      if (.not.ex) then
        open(unit3,file=filename(1:len_trim(filename)), 	 &
  	       form='formatted',				 &
  	       access='SEQUENTIAL', 			  	 &
  	       position='REWIND',				 &
  	       status='unknown')
      endif
!
!      call cloud_stats ( clouds, clouds_old, updrafts, tav, nav, navb, 		&
!    		         nupav, nupavb, nvar, nvarb, nupvar, nupvarb, pdf, 		&
!		         mav, mavb, mvar, mvarb, pdfm, mupav, mupavb, mupvar, mupvarb, 	&
!		         msummax, mupsummax, pdfupm, dav, dvar, 			&
!		         pdfd, dupav, dupvar, pdfupd, pdft, pdfdq, pdfda, pdfdd )
!
!  Outputs
!
!      call cloud_output ( unit1, unit2, unit3, time(t), tav, maxval(clouds%size),	&
!    			  nav, nvar, navb, nvarb, nupav, nupvar, 			&
!			  nupavb, nupvarb, mav, mvar, mavb, mvarb, 			&
!			  mupav, mupvar, mupavb, mupvarb, msummax, mupsummax, 		&
!			  dav, dvar, dupav, dupvar, pdf, pdfm, pdfupm, pdfd, pdfupd, 	&
!			  pdft, pdfdq, pdfda, pdfdd )
!
      close (unit1)
      close (unit2)
      close (unit3)
!
    endif
!
    if (lfirst) lfirst=.false.
    
    if (nt > 0.) exit
!
    endif
!
  enddo
!
!----------------------------------------------------------!
!		     Deallocate and end                    !
!----------------------------------------------------------!
!
  deallocate( mom_all, qt_all, dq_all, buo_all, bef_all, mse_all, det_all, ent_all, liq_all )
  deallocate( cl_id, cl_rad, cl_mfl, cl_mse, cl_buo, cl_ent, cl_det, up_id )
  deallocate( if1, if2 )
  deallocate( nupdrafts )
!
  deallocate( xx, yy, zz, time )
  deallocate( wind, state )
!
!----------------------------------------------------------!
!  
  end program cloud_tracking
