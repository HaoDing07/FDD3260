#include "ctrparam.h"

module surface_var

 !General modules

    USE gridno 
    USE shared_data    
    USE shared_all
    USE averages
    USE gradients
  
    IMPLICIT NONE

    public :: surface_forcing
    
    contains

    !  ===================================================
     subroutine surface_forcing
    !  ===================================================
    !
      call read_forcing_file
      if(isurf == 0) call read_forcing_surface_lhf
      if(isurf == 0) call read_forcing_surface_shf
      call forcing_surface_sst    !
      if(isurf == 0) call forcing_surface_hf

      call read_forcing_surface_ssm
      call forcing_surface_ssm

    return
    end

    !======================================================
    subroutine read_forcing_file
    !=======================================================
      ! USE netcdf
      USE shared_data 
      USE shared_all   
      USE shared_aerosol_new
      USE shared_wind

      USE shared_state
      USE averages
      USE gridno
      USE shared_thermo
      USE gradients
      ! USE mpi
      
      IMPLICIT NONE

      integer :: time_simulation, i, j, m, n

      if (verbose > 0) write(7,*) 'Starting subroutine read_forcing_surface'      
      if (verbose > 0) write(7,*) 'SST forcing file is:',forcing_file 

    ! Read the sst forcing file
      open(233, file=forcing_file, status='old')
      n = 0
    do
        read(233,*,iostat=i)  ! read the data of next line
        if (i /= 0) exit     ! if no data, quit the loop
        n = n + 1            
    end do

    if (.not. allocated(ts_force)) then 
      allocate(ts_force(n))       ! allocate the array
    end if
    rewind(233)              ! let the pointer be at the start of the file

    do i = 1, n
        read(233,*) ts_force(i)   ! store the data in ts_force array
    end do
      close(233)

    if (verbose > 0) write(7,*) 'Reading data and closing forcing file done'

    
    return
    end
     
   !======================================================
    subroutine read_forcing_surface_shf
      !=======================================================
      ! USE netcdf
      USE shared_data 
      USE shared_all   
      USE shared_aerosol_new
      USE shared_wind
      
      USE shared_state
      USE averages
      USE gridno
      USE shared_thermo
      ! USE mpi
      
      IMPLICIT NONE
    
      integer :: time_simulation, i, j, m, n
      if (isurf==0 .and. verbose > 0) write(7,*) 'Starting subroutine read_forcing_surface_shf'            
      if (isurf==0 .and. verbose > 0) write(7,*) 'SHF forcing file is:',forcing_file_shf
      ! Read the shf forcing file
      open(244, file=forcing_file_shf, status='old')
      m = 0
    do
        read(244,*,iostat=i)  ! read the data of next line
        if (i /= 0) exit     ! if no data, quit the loop
        m = m + 1            
    end do

    if (.not. allocated(shf_force)) then 
      allocate(shf_force(m))       ! allocate the array
    end if
    rewind(244)              ! let the pointer be at the start of the file
    
    do i = 1, m
        read(244,*) shf_force(i)   ! store the data in ts_force array
    end do
      close(244)
    
      if (isurf==0 .and. verbose > 0) write(7,*) 'Reading data and closing shf forcing file done'
    
      return
      end
    
   !======================================================
      subroutine read_forcing_surface_lhf
        !=======================================================
      ! USE netcdf
      USE shared_data 
      USE shared_all   
      USE shared_aerosol_new
      USE shared_wind

      USE shared_state
      USE averages
      USE gridno
      USE shared_thermo
      ! USE mpi
        
      IMPLICIT NONE
      
      integer :: time_simulation, i, j, m, n
      if (isurf==0 .and. verbose > 0) write(7,*) 'Starting subroutine read_forcing_surface_lhf'          
      if (isurf==0 .and. verbose > 0) write(7,*) 'LHF forcing file is:',forcing_file_lhf
      ! Read the lhf forcing file
      open(255, file=forcing_file_lhf, status='old')
      j = 0
    do
        read(255,*,iostat=i)  ! read the data of next line
        if (i /= 0) exit     ! if no data, quit the loop
        j = j + 1            
    end do
    
    if (.not. allocated(lhf_force)) then 
      allocate(lhf_force(j))       ! allocate the array
    end if
    
    rewind(255)              ! let the pointer be at the start of the file
    
    do i = 1, j
        read(255,*) lhf_force(i)   ! store the data in ts_force array
    end do
    
      close(255)
    
      if (isurf==0 .and. verbose > 0) write(7,*) 'Reading data and closing lhf forcing file done' 
    
      return
      end
    

    !======================================================
    subroutine forcing_surface_sst
    !=======================================================
      USE shared_state
      USE averages
      USE gridno
      USE shared_thermo
      
      IMPLICIT NONE
      integer :: time_simulation

      if (verbose > 0) write(7,*) 'Starting subroutine forcing_surface_sst'

      time_simulation = 1+(time/60)
    
      if (verbose > 0) write(7,*) 'time_simulation =', time_simulation

      sst = MOD(time,60.)* (ts_force(time_simulation+1) - ts_force(time_simulation))/60. + ts_force(time_simulation)
      
      if (verbose > 0) write(7,*) 'sst =', sst

      if (verbose > 0) write(7,*) 'Finishing subroutine forcing_surface_sst' 

      return
      end

    !======================================================
      subroutine forcing_surface_hf
    !=======================================================
      USE shared_state
      USE averages
      USE gridno
      USE shared_thermo
      
      IMPLICIT NONE
      integer :: time_simulation_hf

      if (isurf==0 .and. verbose > 0) write(7,*) 'Starting subroutine forcing_surface_hf'

      time_simulation_hf = 1+(time/3600)

      if (isurf==0 .and. verbose > 0) write(7,*) 'time_simulation_hf =', time_simulation_hf

      shf = MOD(time,3600.)* (shf_force(time_simulation_hf+1) - shf_force(time_simulation_hf))/3600. + shf_force(time_simulation_hf)
     
      lhf = MOD(time,3600.)* (lhf_force(time_simulation_hf+1) - lhf_force(time_simulation_hf))/3600. + lhf_force(time_simulation_hf)
     
      if (isurf==0 .and. verbose > 0)  write(7,*) 'shf =', shf
      if (isurf==0 .and. verbose > 0)  write(7,*) 'lhf =', lhf
      if (isurf==0 .and. verbose > 0)  write(7,*) 'Finishing subroutine forcing_surface_hf'

      return
      end

  !======================================================
      subroutine read_forcing_surface_ssm
        !=======================================================
          ! USE netcdf
          USE shared_data 
          USE shared_all   
          USE shared_aerosol_new
          USE shared_wind
    
          USE shared_state
          USE averages
          USE gridno
          USE shared_thermo
          USE gradients
          ! USE mpi
          
          IMPLICIT NONE
    
          integer :: time_simulation, i, j, m, n
    
          if (verbose > 0) write(7,*) 'Starting subroutine read_forcing_surface_ssm'      
          if (verbose > 0) write(7,*) 'SST forcing file is:',forcing_file_ssm
    
        ! Read the sst forcing file
          open(266, file=forcing_file_ssm, status='old')
          n = 0
        do
            read(266,*,iostat=i)  ! read the data of next line
            if (i /= 0) exit     ! if no data, quit the loop
            n = n + 1            
        end do
    
        if (.not. allocated(ssm_force)) then 
          allocate(ssm_force(n))       ! allocate the array
        end if
        rewind(266)              ! let the pointer be at the start of the file
    
        do i = 1, n
            read(266,*) ssm_force(i)   ! store the data in ts_force array
        end do
          close(266)
    
        if (verbose > 0) write(7,*) 'Reading data and closing forcing ssm file done'
    
        
        return
      end subroutine read_forcing_surface_ssm

    !======================================================
      subroutine forcing_surface_ssm
        !=======================================================
          USE shared_state
          USE averages
          USE gridno
          USE shared_thermo
          
          IMPLICIT NONE
          integer :: time_simulation
    
          if (verbose > 0) write(7,*) 'Starting subroutine forcing_surface_ssm'
    
          time_simulation = 1+(time/60)
        
          if (verbose > 0) write(7,*) 'time_simulation =', time_simulation
    
          ssm = MOD(time,60.)* (ssm_force(time_simulation+1) - ssm_force(time_simulation))/60. + ssm_force(time_simulation)
          
          if (verbose > 0) write(7,*) 'ssm =', ssm
    
          if (verbose > 0) write(7,*) 'Finishing subroutine forcing_surface_ssm' 
    
          return
          end subroutine forcing_surface_ssm

end module surface_var