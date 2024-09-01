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
!  OUTPUT.F                   
!
!  Purpose:
!	A package of routines used by output.f90			  
!
!  Author
!	Julien Savre, MISU
!
! ================================================================

!	=================================================
	subroutine outx ( name, x, cnt, flag )
!	=================================================

#ifdef SPMD
        USE mpi
#endif
	USE gridno
	USE shared_data
#ifdef SPMD
	USE mpicomm
#endif
	
	IMPLICIT NONE

	real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: x
#ifdef PARALLEL_OUT
	real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: p_x
#else
	real, dimension(nx,ny,nz) :: p_x
#endif
!
    integer  :: flag, unit, i1, i2, j1, j2, nlg, ierr
    real     :: cnt, data_max, data_min, maxall, minall
    character(len=*) :: name
!
#ifdef PARALLEL_OUT
  i1 = 3
  i2 = nx-2
#ifdef MODEL_3D
  j1 = 3
  j2 = ny-2
#else
  j1 = 1
  j2 = 1
#endif
#else
  i1 = it_start
  i2 = it_end
#ifdef MODEL_3D
  j1 = jt_start
  j2 = jt_end
#else
  j1 = 1
  j2 = 1
#endif
#endif
!
  unit = 100
  nlg = nz*(i2-i1+1)*(j2-j1+1)
!
#if ( defined SPMD ) && ( !defined PARALLEL_OUT )
  call collect ( x,  p_x )
#else  
  p_x = x
#endif
!
  p_x = cnt * p_x
!
!  Find and print out min max values
!
#if ( defined SPMD ) && ( !defined PARALLEL_OUT )
        data_max  = maxval (p_x(i1:i2,j1:j2,1:nz))
        data_min  = minval (p_x(i1:i2,j1:j2,1:nz))
!
        CALL MPI_ALLREDUCE (data_max, maxall, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr) 
        CALL MPI_ALLREDUCE (data_min, minall, 1, REALTYPE, MPI_MIN, MPI_COMM_WORLD, ierr) 
#else
        if (mypid == 0) maxall = maxval (p_x(i1:i2,j1:j2,1:nz))
        if (mypid == 0) minall = minval (p_x(i1:i2,j1:j2,1:nz))
#endif
!
       if ( mypid == 0 ) then     
         if (flag == 1)  write(7,9100)trim(name)//' ',maxall,minall
       endif
!
!  Write 3D field
!
#ifdef PARALLEL_OUT
       call open_out (name, unit)
	   call write_to_file_r4(unit,nlg,1.0,p_x(i1:i2,j1:j2,1:nz))		
	   close (unit)        
#else
	   if (mypid==0) then
         call open_out (name, unit)
	     call write_to_file_r4(unit,nlg,1.0,p_x(i1:i2,j1:j2,1:nz))
	     close (unit)        
	   endif
#endif
!
9100	format(11x,"Max & Min ", a6,"      ",f13.4,f16.4)
!
return 
end
!
!	=================================================
	subroutine open_out (vname,fid)
!	=================================================

    USE gridno
    USE shared_data
	
    integer	       :: fid, i
    logical	       :: ex1,ex2
    character(len=*)   :: vname
    character(len=100) :: filename
    character(len=4)   :: cpid
    character(len=5)   :: tcar5
    character(len=6)   :: tcar6
!
!  File name
!
    write(tcar5,'(i5)') floor(time/dt0_i)
    tcar6 = '.' // adjustl(trim(tcar5))
    filename = opdir(1:len_trim(opdir)) // trim(adjustl(vname)) // trim(adjustl(tcar6))
!
!  Inquire file and unit
!
        inquire (file=trim(filename),exist=ex1)
!
	    i = 0
	    inquire (unit=fid,exist=ex2)
	    do while (.not.ex2)
	      fid=fid+i
	      inquire (unit=fid,exist=ex2)
	      i=i+1
	    enddo
!
!  Possibility to write one output file per proc.
!
#ifdef PARALLEL_OUT
    write(cpid,'(a4)') mypid
	open(fid,file=filename(1:len_trim(filename)) // '_' // trim(adjustl(cpid)),   &
     	       form='formatted',		 	&
     	       access='SEQUENTIAL',			&
     	       position=how_access,			&
     	       status='unknown')	
#else
    if (mypid == 0) then
	  open(fid,file=filename(1:len_trim(filename)),   &
     	       form='formatted',		 	&
     	       access='SEQUENTIAL',			&
     	       position=how_access,			&
     	       status='unknown')	    
    endif
#endif

	return
	end subroutine open_out
!
!	=================================================
	subroutine output_ts ( unit_number, ts_tab )	
!	=================================================

	USE gridno
	USE shared_data
!
	integer :: unit_number
	real, dimension(ntsout) :: ts_tab

! ============================================================
!
!  Data
!
  if (trim(ts_out) == 'stratus') then
    if (new_run .and. time == 0.) write(unit_number,100)
    write(unit_number,101)time,ts_tab(1:20)
  else if (trim(ts_out) == 'cumulus') then
    if (new_run .and. time == 0.) write(unit_number,102)
    write(unit_number,103)time,ts_tab(1:32)
  else if (trim(ts_out) == 'deep') then
    if (new_run .and. time == 0.) write(unit_number,104)
    write(unit_number,105)time,ts_tab(1:19)
  else if (trim(ts_out) == 'conservation') then
    if (new_run .and. time == 0.) write(unit_number,106)
    write(unit_number,107)time,ts_tab(1:11)
  endif
!
100  format ("#",3x,"time", 14x,"LWP",                &
         12x,"IWP", 13x,"Zi", 12x,"c_base",           &
         9x,"c_frac", 7x,"surf_prec",                 &
         8x,"max_w",9x,"max_w2",9x,"BIR",	      & 
         10x,"LWP_var", 9x,"Zi_var", 11x,"Zb_var",    &
         12x, "E", 12x, "max_dep", 10x, "SSR", 12x,   &
	 "OLR", 12x, "OLR_CS", 9x, "TOA", 12x,        &
	 "TOA_CS", 9x, "Nc_mean" )
101  format (2x, 21e15.7)
!
102  format ("#",3x,"time", 12x,"max_w",11x,"mean_ct",9x,"mean_cb",9x,"max_cc",     &
             11x,"CF_Z1",10x,"MF_Z1",9x,"CCF_Z1",9x,"CMF_Z1",           &
             10x,"CF_Z2",10x,"MF_Z2",9x,"CCF_Z2",9x,"CMF_Z2",           &
             10x,"CF_Z3",10x,"MF_Z3",9x,"CCF_Z3",9x,"CMF_Z3",           &
             10x,"CWP",11x,"PW",12x,"CRH",11x,"var_PW",                 &
             8x,"prec",10x,"max_prec",7x,"peff",8x,"SHF",11x,"LHF",   	&
             11x,"OLR",11x,"OLR_CS",8x,"TOA",11x,"TOA_CS",              &
             8x,"CAPE",10x,"CIN",11x,"LFC" )
103  format (2x, 33e15.7)
!
104  format ("#",3x,"time", 12x,"max_w",                        &
             11x,"mean_ct",9x,"max_ct",10x,"mean_cb",           &
             9x,"c_frac",10x,"prec",8x,"max_prec",              &
             9x,"O3_trop",8x,"T_trop",9x,"Z_trop",              &
             9x,"Int_O3",9x,"SHF",11x,"LHF",11x,"RADF",         &
             11x,"OLR",11x,"OLR_CS",9x,"TOA",11x,"TOA_CS" )
105  format (2x, 20e15.7)
!
106  format ("#",3x,"time", 12x,"max_w",14x,"LWP",12x,"PW",     &
             13x,"total_dpt",6x,"total_drho",5x,"total_dqt",    &
             6x,"den_tot",8x,"qt_tot",9x,"pt_tot",              &
             8x,"SHF",11x,"LHF")
107  format (2x, 12e15.7)
!	
RETURN
END
!
!	===================================================
	subroutine write_to_prt(unit_number,line,acces)
!	===================================================

        USE shared_data
!
!  	A subroutine to convert data using factor aw_x 
!		then write them as real*4 into a file
!	
	character(len=*)      :: line
	character(len=*)      :: acces
	integer		      :: unit_number

	if (mypid == 0) then
          open(7,file=prtname(1:len_trim(prtname)), access='SEQUENTIAL', position=trim(acces), status='unknown')

          write(7,*) trim(line)

          close(7)
	endif

100	format(6e15.7)
	
	return
	end
!
!	===================================================
	subroutine write_to_file_r4(unit_number,nl,aw_x,x)
!	===================================================
!
!  	A subroutine to convert data using factor aw_x 
!		then write them as real*4 into a file
!	
	real(8),dimension(nl) :: spep
	real,   dimension(nl) :: x
	integer		      :: nl, unit_number
	real		      :: aw_x

	spep = x
	
	if (unit_number < 1000) rewind(unit_number)

	write(unit_number,100)spep

100	format(6e15.7)
	
	return
	end
!
!	=======================================================
	subroutine repositionts (unit)
!	=======================================================
	
	USE shared_data

	integer		 :: unit, io
	character(len=1) :: car1
	
	!
	!  Test end of file
	!
	read(unit,'(a1)',iostat=io) car1

    	do while (io == 0)
	  read(unit,'(a1)',iostat=io) car1
          how_access = 'append'
    	enddo
	backspace (unit)
	
	return
	end subroutine repositionts	
!	
!	=======================================================
	subroutine reposition (unit, nd)
!	=======================================================
	
	USE shared_data

	integer		 :: i, j, io, unit, flag, nl, nd
	real		 :: timep
	
	read(unit,'(a1)',iostat=io)
	rewind(unit)
	
    	if (io == 0 .and. .not.new_run) then
    
	nl = ceiling(real(nz)/6.)
	timep = time+dt
	flag = 0	
	do while (time - timep < 0.)
	  flag = flag + 1
          do i = 1, nd
	    do j = 1, nl
	      backspace (unit)
	    enddo
	    if (nd /= 1) backspace (unit)
	  enddo	  
	  backspace (unit)	  
	  read(unit,*) timep
	  backspace (unit)
	  backspace (unit)
	enddo
	read(unit,*)
	read(unit,*)

	do i = 1, nd
	  if (nd /= 1) read(unit,*)
	  do j = 1, nl
	    read(unit,*)
	  enddo
    	enddo
    
    	endif

	return
	end subroutine reposition	
