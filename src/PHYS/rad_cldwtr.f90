!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module cldwtr

  use gridno, only : nhab, datadir, nz
  use shared_data, only : ihab, nmode, mypid, pi
  use radia, only : norig, nv, mb, mbs, mbir
  implicit none
  integer, save :: nsizes
  logical, save :: Initialized = .False.
  LOGICAL, SAVE :: aeroInitialized = .FALSE.  !UCLALES-SALSA

  real, allocatable    :: re(:), fl(:), bz(:,:), wz(:,:), gz(:,:)
  real, allocatable    :: brn(:), wrnf(:), grn(:)
  real, allocatable    :: ap(:,:), bp(:,:), cps(:,:,:), dps(:,:), cpir(:,:)
  real, allocatable    :: aip(:,:,:), bip(:,:,:), cipgf(:,:,:)

  ! Juha (UCLALES-SALSA):  
  ! Lookup table variables for aerosol optical properties for radiation calculations
  ! Real and imaginary parts of refractive indices, size parameter, extinction crossection, asymmetry parameter and omega???
  REAL, ALLOCATABLE, TARGET :: aer_nre_LW(:), aer_nim_LW(:), aer_alpha_LW(:),   &
                       aer_sigma_LW(:,:,:), aer_asym_LW(:,:,:), aer_omega_LW(:,:,:)
  REAL, ALLOCATABLE, TARGET :: aer_nre_SW(:), aer_nim_SW(:), aer_alpha_SW(:),   &
                       aer_sigma_SW(:,:,:), aer_asym_SW(:,:,:), aer_omega_SW(:,:,:)
  ! /Juha


contains
  !
  !---------------------------------------------------------------------------
  ! Surbourine cloud_init initialize data arrays for the cloud model, 
  ! checking for consistency between band structure of cloud model and CKD
  !
  subroutine init_cldwtr

    use ckd, only : band, center
    integer, parameter  :: nrec = 21600

    real, dimension(mb) :: cntrs

    integer             :: ib, i, j, k, nbands, n1, n2
    character (len=12)  :: frmt

    open ( unit = 41, file = datadir(1:len_trim(datadir))//'/cldwtr.dat', status = 'old', recl=nrec )
    open ( unit = 42, file = datadir(1:len_trim(datadir))//'/cirrwtr.dat', status = 'old', recl=nrec )
    open ( unit = 43, file = datadir(1:len_trim(datadir))//'/rnwtr.dat', status = 'old', recl=nrec )
    open ( unit = 44, file = datadir(1:len_trim(datadir))//'/icewtr.dat', status = 'old', recl=nrec )

    read (41,'(2I3)') nsizes, nbands
    if (nbands /= mb .or. nsizes*nbands*15 > nrec) &
         stop 'TERMINATING: incompatible cldwtr.dat file'

    allocate (re(nsizes),fl(nsizes),bz(nsizes,mb),wz(nsizes,mb),gz(nsizes,mb))
    write(frmt,'(A1,I2.2,A8)') '(',mb,'E15.7)    '
    read (41,frmt) (cntrs(i), i=1,mb)
    do i=1,mb
       if (spacing(1.) < abs(cntrs(i) - center(band(i))) ) &
            stop 'TERMINATING: cloud properties not matched to band structure'
    end do

    write(frmt,'(A1,I2.2,A9)') '(',nsizes,'E15.7)   '
    read (41,frmt) (re(i), i=1,nsizes)
    read (41,frmt) (fl(i), i=1,nsizes)

    write(frmt,'(A1,I4.4,A7)') '(',nsizes*mb,'E15.7) '
    read (41,frmt) ((bz(i,ib), i=1,nsizes), ib=1,mb)
    read (41,frmt) ((wz(i,ib), i=1,nsizes), ib=1,mb)
    read (41,frmt) ((gz(i,ib), i=1,nsizes), ib=1,mb)
    close (41)

    read (42,*) n1, n2
    if (mbs /= n1 .or. mbir /= n2) &
         stop 'TERMINATING: incompatible cirrwtr.dat file'
    allocate (ap(3,mb), bp(4,mb), cps(4,4,mbs), dps(4,mbs), cpir(4,mbir))        
	         
    read(42,*) ((ap(i,ib), i=1,3), ib=1,mb)
    read(42,*) ((bp(i,ib), i=1,4), ib=1,mb)
    
    read(42,*) (((cps(i,j,ib), i=1,4), j=1,4), ib=1,mbs)
    read(42,*) ((dps(i,ib), i=1,4), ib=1,mbs)
    read(42,*) ((cpir(i,ib), i=1,4), ib=1,mbir)
    close (42)

    read (43,*) nbands
    if (nbands /= mb) &
         stop 'TERMINATING: incompatible cirrwtr.dat file'
    allocate (brn(mb), wrnf(mb), grn(mb))

    read(43,*) (brn(ib), ib=1,mb)
    read(43,*) (wrnf(ib), ib=1,mb)
    read(43,*) (grn(ib), ib=1,mb)
    close (43)

    read (44,*) n1, n2
    if (mbs /= n1 .or. nhab /= n2) &
         stop 'TERMINATING: incompatible icewtr.dat file'
    allocate (aip(nhab,4,mbs), bip(nhab,4,mbs), cipgf(nhab,6,mbs))        
     
    do k = 1, nhab
      read(44,*) ((aip(k,i,ib), i=1,4), ib=1,mbs)
      read(44,*) ((bip(k,i,ib), i=1,4), ib=1,mbs)
      read(44,*) ((cipgf(k,i,ib), i=1,6), ib=1,mbs)
    enddo
    close (44)
        
    if (minval((/bz,wz,gz/)) < 0.) &
         stop 'TERMINATING: cloud properties out of bounds'

    Initialized = .True.

  end subroutine init_cldwtr



  ! Taken from UCLALES-SALSA
  !---------------------------------------------------------------
  ! ReadAeroRadInput
  !
  SUBROUTINE init_aerorad
   IMPLICIT NONE

   ! Shortwave tables
   CALL init_aerorad_lookuptables(datadir(1:len_trim(datadir))//'/lut_uclales_salsa_sw.nc', & 
                                  aer_nre_SW, aer_nim_SW, aer_alpha_SW,  &
                                  aer_sigma_SW, aer_asym_SW, aer_omega_SW          )

   
   ! Longwave tables
   CALL init_aerorad_lookuptables(datadir(1:len_trim(datadir))//'/lut_uclales_salsa_lw.nc', &
                                  aer_nre_LW, aer_nim_LW, aer_alpha_LW,  &
                                  aer_sigma_LW, aer_asym_LW, aer_omega_LW          )


   aeroInitialized = .TRUE.



 END SUBROUTINE init_aerorad

 !Taken from UCLALES-SALSA
 SUBROUTINE init_aerorad_lookuptables(filename, zaer_nre, zaer_nim, zaer_alpha, &
                                      zaer_sigma, zaer_asym, zaer_omega         )
   USE mo_aerorad_lut, ONLY : init_aerorad_lut
   IMPLICIT NONE

   CHARACTER(len=*), INTENT(in) :: filename
   REAL, ALLOCATABLE, INTENT(out) :: zaer_nre(:), zaer_nim(:), zaer_alpha(:)
   REAL, ALLOCATABLE, INTENT(out) :: zaer_sigma(:,:,:), zaer_asym(:,:,:), zaer_omega(:,:,:)

   INTEGER :: Nalpha, Nim, Nre

   CALL init_aerorad_lut(filename, Nalpha, Nim, Nre,        &
   zaer_nre, zaer_nim, zaer_alpha,    &
   zaer_sigma, zaer_asym, zaer_omega  )

 END SUBROUTINE init_aerorad_lookuptables



  ! -----------------------------------------------------------------------
  ! Subroutine cloud_water:  calculates the optical depth (tw), single 
  ! scattering albedo (ww), and phase function (www(4)) given the cloud 
  ! water [g/m^3] and effective radius [microns] by interpolating based on
  ! known optical properties at predefined sizes  
  !
  subroutine cloud_water ( ib, pre, pcw, dz, tw, ww, www )

    implicit none

    integer, intent (in) :: ib
    real, dimension (nv), intent (in) :: pre, pcw, dz
    real, intent (out) :: tw(nv), ww(nv), www(nv,4)

    integer :: k, j, j0, j1
    real    :: gg, wght, cwmks
!
    if (.not.Initialized) stop 'TERMINATING: Cloud not Initialized'
!
    www = 0.0
    tw  = 0.0
    ww  = 0.0
    do k = 1, nv
       cwmks = pcw(k)*1.e-3
       if ( cwmks .ge. 1.e-8) then
          j = 0
          do while (pre(k) > re(j+1))
             j = j + 1
	     if (j == nsizes) exit
          end do
          if (j >= 1 .and. j < nsizes) then
             j1 = j+1
             wght = (pre(k)-re(j))/(re(j1)-re(j))
             tw(k) = dz(k) * cwmks * ( bz(j,ib) / fl(j) +   &
                  ( bz(j1,ib) / fl(j1) - bz(j,ib) / fl(j) ) /    &
                  ( 1.0 / re(j1) - 1.0 / re(j) ) * ( 1.0 / pre(k) &
                  - 1.0 / re(j) ) )
             ww(k) = wz(j,ib) + (wz(j1,ib) - wz(j,ib) ) * wght
             gg    = gz(j,ib) + (gz(j1,ib) - gz(j,ib) ) * wght
          else
             j0 = max(j,1)
             tw(k) = dz(k) * cwmks * (bz(j0,ib)/fl(j0))
             ww(k) = wz(j0,ib) 
             gg    = gz(j0,ib)          
          end if
          www(k,1) = 3.0 * gg
          do j=2,4
             wght = real(2*j+1)/real(2*j-1)
             www(k,j) = www(k,j-1) * gg * wght
          end do
       end if
    end do

    return
  end subroutine cloud_water


  ! -----------------------------------------------------------------------
  ! Subroutine rain_water:  calculates the optical depth (tw), single 
  ! scattering albedo (ww), and phase function (www(4)) given the rain 
  ! water [g/m^3] by interpolating based on
  ! known optical properties at predefined sizes  
  !
  subroutine rain_water ( ib, prw, dz, tw, ww, www )

    implicit none

    integer, intent (in) :: ib
    real, dimension (nv), intent (in) :: prw, dz
    real, intent (out) :: tw(nv), ww(nv), www(nv,4)

    integer :: j, k
    real    :: wght, cwmks, rwc=0.5
!
!  Warning: computation of rain optical properties is only performed 
!	    inside the LES domain.
!
!
    www = 0.0
    tw  = 0.0
    ww  = 0.0
    do k = norig+1, nv
       cwmks = prw(k)*1.0e-3
       if ( cwmks .ge. 1.0e-8 ) then
    	 tw(k) = dz(k)*cwmks * brn(ib) / rwc
    	 ww(k) = wrnf(ib)
!
         www(k,1) = 3.0 * grn(ib)
         do j=2,4
            wght = real(2*j+1)/real(2*j-1)
            www(k,j) = www(k,j-1) * grn(ib) * wght
         end do
       endif
    enddo

    return
  end subroutine rain_water
  
  ! -----------------------------------------------------------------------
  ! Subroutine ice_water:  calculates the optical depth (tw), single 
  ! scattering albedo (ww), and phase function (www(4)) given the ice 
  ! water [g/m^3] and effective radius [microns] by interpolating based on
  ! known optical properties at predefined sizes  
  !
  subroutine ice_water ( ib, pde, piwc, dz, tw, ww, www )

    implicit none

    integer, intent (in) :: ib
    real, dimension (nv), intent (in) :: pde, piwc, dz
    real, intent (out) :: tw(nv), ww(nv), www(nv,4)
    
    integer  :: j, k
    real  :: wf, gg, fw1, fw2, fw3, wght, cwmks

    if (.not.Initialized) stop 'TERMINATING: Cloud not Initialized'
!
!  Warning: computation of ice optical properties above the LES domain, 
!	    that is for cirrus clouds only (Parameterization from Fu & Liou).
!   
    www = 0.0
    tw  = 0.0
    ww  = 0.0 
    do k = 1, norig
       cwmks = piwc(k)*1.e-3		    !convert to kg/m3 (input is in g/m3)
       if ( cwmks .ge. 1.e-8 ) then
    	 fw1 = pde(k)
    	 fw2 = fw1 * pde(k)
    	 fw3 = fw2 * pde(k)
    	 tw(k) = 1000.*dz(k)*cwmks * (ap(1,ib) + ap(2,ib) / fw1 + ap(3,ib) / fw2)
    	 ww(k) = 1.0 - bp(1,ib) + bp(2,ib) * fw1 + bp(3,ib) * fw2 + bp(4,ib) * fw3
    	 if ( ib .le. mbs ) then
    	   gg = dps(1,ib) + dps(2,ib) * fw1 + dps(3,ib) * fw2 + dps(4,ib) * fw3
    	   do j=1,4
    	      wght = real(2*j+1)
    	      wf = cps(1,j,ib) + cps(2,j,ib)*fw1 + cps(3,j,ib)*fw2 + cps(4,j,ib)*fw3
    	      www(k,j) = ( 1.0 - gg )*wf + wght*gg
    	   end do
    	 else
    	   gg = cpir(1,ib-mbs) + cpir(2,ib-mbs)*fw1 + cpir(3,ib-mbs)*fw2 + cpir(4,ib-mbs)*fw3
    	   www(k,1) = 3.0 * gg
    	   do j=2,4
    	      wght = real(2*j+1)/real(2*j-1)
    	      www(k,j) = www(k,j-1) * gg * wght
    	   end do
    	 endif
       endif
    end do
!
!  Warning: Now computation of ice optical properties inside the LES domain, for
!	    low-level ice clouds (Parameterization for shortwave only from Key et al., 2002)
! 
!    do k = norig+1, nv
!    	cwmks = piwc(k)*1.e-3		     !convert to kg/m3 (input is in g/m3)
!    	if ( cwmks .ge. 1.e-8 .and. ib .le. mbs) then
!    	  fw1 = pde(k)
!    	  fw2 = fw1 * pde(k)
!    	  fw3 = fw2 * pde(k)
!    	  tw(k) = dz(k)*cwmks * (aip(ihab,1,ib) + aip(ihab,2,ib) / fw1 + aip(ihab,3,ib) / fw2 + aip(ihab,4,ib) / fw3)
!    	  ww(k) = bip(ihab,1,ib) + bip(ihab,2,ib) * fw1 + bip(ihab,3,ib) * fw2 + bip(ihab,4,ib) * fw3
!    	  gg	= cipgf(ihab,1,ib) + cipgf(ihab,2,ib) * fw1 + cipgf(ihab,3,ib) * fw2 + cipgf(ihab,4,ib) * fw3

!    	   www(k,1) = 3.0 * gg
!    	   do j=2,4
!    	     wght = real(2*j+1)/real(2*j-1)
!    	     www(k,j) = www(k,j-1) * gg * wght
!    	   end do
!    	 endif
!    end do
!	
    return
  end subroutine ice_water
 
!
!Taken from UCLA-LES Salsa with some modifications 
! Calculates the optical depth (taer), single scattering albedo (waer) and phase function (wwaer(4)) for given
  ! binned aerosol mass and number concentration arrays using lookup tables for optical properties.
  SUBROUTINE aero_rad(ib, volc,voltot, aero1d_num, dz, taer, waer, wwaer)
   USE ckd, ONLY : band, center, IsSolar,llimit,rlimit 
   USE util  !Here we find functions getMAssIndex and closest, originally from module "util" in SALSA.  
   USE mo_salsa_optical_properties, ONLY : aerRefrIBands_SW, aerRefrIBands_LW,  &
                                           riReSW, riImSW, riReLW, riImLW
   USE shared_aerosol_new, only : nelem, aero
   
   IMPLICIT NONE

   INTEGER, INTENT(in) :: ib            !number of bands
   REAL, PARAMETER   :: nlim=1.         ! Number conc. limit (#/kg) for aerosol and cloud droplets  (from mo_submctl.f90)  
   REAL, PARAMETER   :: pi6=0.5235988   !pi/6 taken from mo_submctl.f90 in SALSA
   real, dimension (nv), intent (in) :: dz
   real, INTENT(out) :: taer(nv), waer(nv), wwaer(nv,4)   ! optical depth, single scattering albedo, phase function
   
   REAL :: lambda_r   ! Center wavenumber for current band 1/cm
   Real :: lambda_left, lambda_right, lambda_center
   
   real, dimension (nmode,nelem,nv) :: volc        ! Particle volume concentrations for each mode, specie and vertical level
   real, dimension (nmode,nv)       :: voltot, aero1d_num      ! Total particle volume for each mode and vertical level                          


   ! Refractive index for each chemical (closest to current band from tables in submctl)
   REAL              :: refrRe_all(nelem), refrIm_all(nelem)
   REAL              :: volmean_refrRe, volmean_refrIm ! Volume mean refractive index for single bin
   !REAL              :: naerotot                       ! Total number of aerosol particles

   ! Size parameter for given wavelength and size bin
   REAL            :: sizeparam, sizeparam2

   ! Binned optical properties; These will be integrated and normalized for final results
   REAL            :: taer_mode(nv,nmode), waer_mode(nv,nmode), wwaer_mode(nv,nmode,4)
   
   ! Lookup table: Indices in refractive index vectors and sizeparameter vector
   INTEGER         :: i_re, i_im, i_alpha

   INTEGER :: refi_ind  ! index for the vector with refractive indices for each wavelength

   ! Bunch of other idices (for loops)
   INTEGER :: ss,kk, istr,iend, k
   INTEGER :: bb

   REAL :: TH = 1.e-30  !Some kind of threshold...?

   REAL, POINTER :: aer_nre(:) => NULL(), aer_nim(:) => NULL(),          &
                    aer_alpha(:) => NULL(), aer_sigma(:,:,:) => NULL(),  &
                    aer_asym(:,:,:) => NULL(), aer_omega(:,:,:) => NULL()

   !IF (.NOT. aeroInitialized) STOP "TERMINATING: Aerosol for radiation not initialized"
   !print *,"Starting aero_rad"

   taer = 0.
   waer = 0.
   wwaer = 0.
   taer_mode = 0.
   waer_mode = 0.
   wwaer_mode = 0.

   lambda_r= center(band(ib))
   IF (lambda_r > 280.) then
      lambda_left= 1./llimit(band(ib))
      lambda_right=1./rlimit(band(ib))
      lambda_center=0.5*(lambda_left+lambda_right)

      !print *,"center band",lambda_r,"lambda_center",lambda_center,"lambda_right",lambda_right


      !lambda_r = center(band(ib))
      !IF (1./lambda_r > aerRefrIbands_SW(1)) THEN
      IF (lambda_center > Maxval(aerRefrIbands_SW)) THEN

         ! Get the refractive indices from the LW tables for the current band
         !refi_ind = closest(aerRefrIbands_LW,1./lambda_r)
         refi_ind = closest(aerRefrIbands_LW,lambda_center)   

         refrRe_all(:) = riReLW(:,refi_ind)
         refrIm_all(:) = riImLW(:,refi_ind)
         aer_nre => aer_nre_LW(:)
         aer_nim => aer_nim_LW(:)
         aer_alpha => aer_alpha_LW(:)
         aer_sigma => aer_sigma_LW(:,:,:)
         aer_asym => aer_asym_LW(:,:,:)
         aer_omega => aer_omega_LW(:,:,:)

      ELSE
         ! Get the refractive indices from the SW tables fr the current band
         !refi_ind = closest(aerRefrIbands_SW,1./lambda_r)
         refi_ind = closest(aerRefrIbands_SW,lambda_center)

         refrRe_all(:) = riReSW(:,refi_ind)
         refrIm_all(:) = riImSW(:,refi_ind)
         aer_nre => aer_nre_SW(:)
         aer_nim => aer_nim_SW(:)
         aer_alpha => aer_alpha_SW(:)
         aer_sigma => aer_sigma_SW(:,:,:)
         aer_asym => aer_asym_SW(:,:,:)
         aer_omega => aer_omega_SW(:,:,:)      

      END IF

   Else

      ! Get the refractive indices from the LW tables fr the current band
      refi_ind = closest(aerRefrIbands_LW,1./lambda_r)
      !refi_ind = closest(aerRefrIbands_LW,lambda_center)   #in this case lambda_right= infinite 

      refrRe_all(:) = riReLW(:,refi_ind)
      refrIm_all(:) = riImLW(:,refi_ind)
      aer_nre => aer_nre_LW(:)
      aer_nim => aer_nim_LW(:)
      aer_alpha => aer_alpha_LW(:)
      aer_sigma => aer_sigma_LW(:,:,:)
      aer_asym => aer_asym_LW(:,:,:)
      aer_omega => aer_omega_LW(:,:,:)
      
   End if


   ! Loop over levels
   do kk =1, nv

      ! loop over modes
      DO bb = 1,nmode 

         ! Check if empty bin
         IF (voltot(bb,kk) < 1.e-30) CYCLE
         
         ! Volume mean refractive indices in current mode
         volmean_refrRe = SUM(volc(bb,1:nelem,kk)*refrRe_all(1:nelem))/voltot(bb,kk)
         volmean_refrIm = SUM(volc(bb,1:nelem,kk)*refrIm_all(1:nelem))/voltot(bb,kk)
            
         ! size parameter in current mode  (uses the volume mean diameter)
         !sizeparam = pi*1.e2*lambda_r*(((voltot(bb,kk)/aero1d_num(bb,kk))/pi6)**(1./3.))

         !Computes the size parameter using the wavelength at the center of the band. 
         !If center wavenumber < 280 then we use the center wave number to compute the size paremeter 
         ! size parameter in current mode  (uses the volume median diameter) 
         IF (lambda_r < 280.) then
            sizeparam = pi*1.e2*lambda_r*(((voltot(bb,kk)/aero1d_num(bb,kk))/pi6)**(1./3.))*exp(3/2*(log(aero(bb)%size%sigma)*log(aero(bb)%size%sigma)))
         ELSE             
            sizeparam = pi*1.e2*(1/lambda_center)*(((voltot(bb,kk)/aero1d_num(bb,kk))/pi6)**(1./3.))*exp(3/2*(log(aero(bb)%size%sigma)*log(aero(bb)%size%sigma)))
         Endif

         ! Corresponding lookup table indices
         i_re = closest(aer_nre,volmean_refrRe)
         i_im = closest(aer_nim,volmean_refrIm)
         i_alpha = closest(aer_alpha,sizeparam)


         ! Binned optical properties
         ! Optical depth  
         IF (lambda_r < 280.) then           
            taer_mode(kk,bb) = 1.e-6*aero1d_num(bb,kk) * dz(kk)*1.e2 * aer_sigma(i_re,i_im,i_alpha) * (1./lambda_r)**2
         else
            taer_mode(kk,bb) = 1.e-6*aero1d_num(bb,kk) * dz(kk)*1.e2 * aer_sigma(i_re,i_im,i_alpha) * (lambda_center)**2  
            
            !write(7,*)'lamda_center=',lambda_center,'aer_omega=',aer_omega(i_re,i_im,i_alpha)   
         endif

         ! Single scattering albedo
         waer_mode(kk,bb) = taer_mode(kk,bb) * aer_omega(i_re,i_im,i_alpha)

         ! Phase function moments
         wwaer_mode(kk,bb,1) = waer_mode(kk,bb) * 3.*aer_asym(i_re,i_im,i_alpha) 

         wwaer_mode(kk,bb,2) = waer_mode(kk,bb) * 5.*aer_asym(i_re,i_im,i_alpha)**2

         wwaer_mode(kk,bb,3) = waer_mode(kk,bb) * 7.*aer_asym(i_re,i_im,i_alpha)**3

         wwaer_mode(kk,bb,4) = waer_mode(kk,bb) * 9.*aer_asym(i_re,i_im,i_alpha)**4

      END DO

      ! Integrate and normalize
      taer(kk) = SUM(taer_mode(kk,1:nmode))
 
      IF (taer(kk) < TH) THEN
         ! Avoid normalization by zero
         taer(kk) = 0.
         waer(kk) = 0.
         wwaer(kk,1:4) = 0.
      ELSE
         waer(kk) = SUM(waer_mode(kk,1:nmode))/taer(kk)
         wwaer(kk,1:4) = SUM(wwaer_mode(kk,1:nmode,1:4),DIM=1)/waer(kk)
      END IF

   END DO

   !print *, "==Total taer==="
   !print *, taer
   !print *, "==============="

   aer_nre => NULL()
   aer_nim => NULL()
   aer_alpha => NULL()
   aer_sigma => NULL()
   aer_asym => NULL()
   aer_omega => NULL()

   !print *,"Ending aero_rad"

 END SUBROUTINE aero_rad


  ! ---------------------------------------------------------------------------
  ! linear interpolation between two points, returns indicies of the 
  ! interpolation points and weights
  !
  subroutine interpolate(x,ny,y,i1,i2,alpha)

    integer, intent (in) :: ny
    real, intent (in)    :: x, y(ny)

    integer, intent (out) :: i1, i2
    real, intent (out)    :: alpha

    if (y(1) < y(2)) stop 'TERMINATING: band centers increasing'

    i2 = 1
    do while (x < y(i2) .and. i2 < ny)
       i2 = i2+1
    end do
    i1 = max(1,i2-1)
    alpha = 1.

    if(i2.ne.i1) alpha = (x-y(i1))/(y(i2)-y(i1))
    if (alpha <0 .or. alpha >1) print 600, x, y(1), y(ny), alpha

    return

600 format(/'CLOUD_INIT WARNING:  Extrapolating because data out of range', &
         /1x,'x = ',F8.1,', ymax = ',F7.0,', ymin =',F7.0,', alpha = ',F6.3)
  end subroutine interpolate

end module cldwtr
