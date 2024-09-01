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
!  SHAREDATA.F:                   
!
!  Purpose:
!      Module holds shared data, a replacement of old INCLUDE/com.h.                    
!
!  Author
!      Chien Wang
!      MIT Joint Program on Science and Policy of Global Change
!            
! ================================================================

      module shared_data

      use gridno
      use, intrinsic :: iso_c_binding 
            
      SAVE
      
  	  !
	  ! --- parameters
	  !
      integer, parameter :: nract    = 100      ! No. of gaseous reactions
      integer, parameter :: nvaria   = 015      ! No. of gaseous variables
      integer, parameter :: nequil   = 23       ! No. of equilibrium reactions
      integer, parameter :: naqreact = 32       ! No. of aqueous reactions
      integer, parameter :: naqvaria = 33       ! No. of variables in aqueouse
      integer            :: nhydro              ! No. of hydrometeors
      integer		 :: nmode0, nmode	! No. of aerosol modes
      integer		 :: ihab		! Index selected ice habit
      integer            :: nlsod               ! No. of elements in the LSODE system
                  
      real, parameter :: Avogadro = 6.02217e+23
      real, parameter :: f2omega  = 1.4584e-4     ! f = 2*omega in Coriolis 
      real, parameter :: rerg = 8.314    ! (Jmol-1K-1) universal gas constant
      real, parameter :: densw = 0.997e3 ! density of water

      !
      ! --- molecular weight: (g mol-1)
      !
      real, parameter :: awN2    = 28.0134, awO2   = 31.9988
      real, parameter :: awH2O   = 18.0152                 
      real, parameter :: awO     = 15.9994, awO3   = 47.9982
      real, parameter :: awCO    = 28.0104, awCO2  = 44.0098
      real, parameter :: awH     = 01.0079, awHO   = 17.0073
      real, parameter :: awHO2   = 33.0067, awH2O2 = 34.0146
      real, parameter :: awNO    = 30.0061, awNO2  = 46.0055
      real, parameter :: awNO3   = 62.0049, awN2O5 =108.0104
      real, parameter :: awHNO3  = 63.0128
      real, parameter :: awCH4   = 16.0426, awCH3  = 15.0347
      real, parameter :: awCHO   = 29.0183, awCH2O = 30.0262
      real, parameter :: awCH3O  = 31.0341, awCH3O2= 47.0335
      real, parameter :: awCH3O2H= 48.0414
      real, parameter :: awSO2   = 64.0588, awHOSO2= 81.0661
      real, parameter :: awSO3   = 80.0582, awH2SO4= 98.0734
      real, parameter :: awF11   =137.3675, awF12  =120.9054
      real, parameter :: awN2O   = 44.0000, awNH4  = 18.0545 
      real, parameter :: awDMS   = 62.0820, awso4 = 96.0468
      real, parameter :: awC5H8  = 68.1185, awNH3 = 17.0466
      
      !
      ! --- sticking coefficient:
      !
      real, parameter :: alphaCO2   = 0.01
      real, parameter :: alphaH2O2  = 0.18
      real, parameter :: alphaO3    = 2.0e-3
      real, parameter :: alphaHNO3  = 0.11
      real, parameter :: alphaSO2   = 0.13
      real, parameter :: alphaH2SO4 = 0.98
      real, parameter :: alphaNH3   = 0.18
      real, parameter :: alphaEST   = 0.1

      !
      ! --- retention coefficient
      !
      real, parameter :: retent_o3       = 1.0
      real, parameter :: retent_civ    = 0.0       ! to keep mass conservation 
      real, parameter :: retent_h2o2   = 1.0       ! 0.05 Snider & Huang (2000)
      real, parameter :: retent_xnv    = 1.0
      real, parameter :: retent_ch2o   = 1.0
      real, parameter :: retent_ch3o2h = 1.0
      real, parameter :: retent_siv    = 1.0       ! 0.62 Iribarne et al.
      real, parameter :: retent_svi    = 1.0
      real, parameter :: retent_nh4    = 1.0       

      real, parameter :: escape_o3     = 1.0 - retent_o3
      real, parameter :: escape_civ    = 1.0 - retent_civ
      real, parameter :: escape_h2o2   = 1.0 - retent_h2o2
      real, parameter :: escape_xnv    = 1.0 - retent_xnv
      real, parameter :: escape_ch2o   = 1.0 - retent_ch2o
      real, parameter :: escape_ch3o2h = 1.0 - retent_ch3o2h
      real, parameter :: escape_siv    = 1.0 - retent_siv
      real, parameter :: escape_svi    = 1.0 - retent_svi
      real, parameter :: escape_nh4    = 1.0 - retent_nh4

      !
      ! --- heterogeneous reaction possibilities (e.g., Molina et al, JPC, 1996)
      !
      real, parameter :: gamma_o3      = 1.e-6         ! Dlugokencky & Ravishankara (1992)
      real, parameter :: gamma_h2o2      = 1.e-4    ! Myhre & Nielsen (1998), 60%svi
      real, parameter :: gamma_xnv      = 0.2         ! Abbatt (1997)
      real, parameter :: gamma_ch2o      = 0.01     ! Tolbert et al (1996) 67%svi
      real, parameter :: gamma_ch3o2h = 9.e-3    ! Magi et al (1997) droplet
      real, parameter :: gamma_siv      = 8.e-6    ! Clegg and Abbatt (2001)
      real, parameter :: gamma_svi      = 0.43     ! Poschl et al (1998) 303Ksvi
      real, parameter :: gamma_nh4      = 1.e-4    ! Same as for H2O2 

!
! --- Min and Max values for hydrometeors
!
      real, parameter :: qc_min  = 1.e-15, xnc_min = 1.e-1, lc_min = 1.e+9, lc_max = 1.e+14
      real, parameter :: qi_min  = 1.e-15, xni_min = 1.e-1, li_min = 1000., li_max = 1.e+9
      real, parameter :: qr_min  = 1.e-12, xnr_min = 1.e-3, lr_min = 100., lr_max = 1.e+6
      real, parameter :: qg_min  = 1.e-12, xng_min = 1.e-3, lg_min = 10., lg_max = 1.e+6
      real, parameter :: qs_min  = 1.e-12, xns_min = 1.e-3, ls_min = 100., ls_max = 1.e+6
      real, parameter :: qh_min  = 1.e-12, xnh_min = 1.e-3, lh_min = 1., lh_max = 1.e+4
      real, parameter :: ccn_min = 1.e-2,  xin_min = 1.e-2, qin_min = 1.e-18
      real, parameter :: qaero_min = 1.e-18, naero_min = 1.e-2 

!
! ---  Physical constants
!
        real, parameter :: pi = 3.14159265
	real, parameter :: g = 9.81		 
	real, parameter :: pref = 100000.		! Reference pressure
	real, parameter :: kb = 1.38065e-23		! Boltzmann constant
	real, parameter :: hp = 6.626068e-34		! Planck
	real, parameter :: kbhp = 2.08366e10		! = kb/hp
	real, parameter :: Na = 6.022e23   		! Avogadro
	real, parameter :: R  = 8.314      		! Universal gas constant
	real, parameter :: Ra = 287.			! Constant for air
	real, parameter :: Rw = 461.5			! Constant for water vapor
	real, parameter :: nmol = 5.85e18		! Number of water molecule at interface
	real, parameter :: mwm = 2.99e-26		! Mass of a water molecule
	real, parameter :: Mw = 0.018014		! Molar weight of water
	real, parameter :: Ma = 0.0289			! Molar weight of air
!	real, parameter :: c1s = 1.e19			! Concentration of vapor molecules at surface
	real, parameter :: nuv = 6.e12			! Vibration frequency of an adsorbed water molecule
	real, parameter :: cp  = 1004.			! Heat capacity
	real, parameter :: cv  = 717.			! Heat capacity
	real, parameter :: rhow = 1000.			! Liquid water density
	real, parameter :: kappas = 0.71		! Kappa value for (H2SO4)
	real, parameter :: rhoi = 920.			! Ice density
	real, parameter :: Kt = 2.4e-2			! Heat conductivity of air
	real, parameter :: Kl = 0.57			! Heat conductivity of liquid water
        real, parameter :: beta = 0.5			! Exponent from KC bulk aerosol parameterization
        real, dimension(3), parameter :: denp = (/500.,300.,2800./)     ! Particle density
	real, dimension(3), parameter :: Ktp = (/0.15,0.15,0.9/)	! Particle thermal conductivity
!
! --- Output codes
!
        logical  :: out_u, out_v, out_w
	logical  :: out_p, out_es, out_rho
	logical  :: out_t, out_ptv, out_pt, out_mse
	logical  :: out_qv, out_qt, out_vp
	logical  :: out_qc, out_nc, out_dc
	logical  :: out_qr, out_nr, out_dr
	logical  :: out_qi, out_ni, out_di
	logical  :: out_qg, out_ng, out_dg, out_wg
	logical  :: out_qs, out_ns, out_ds, out_ws
	logical  :: out_qh, out_nh, out_dh, out_wh
	logical  :: out_var, out_sca, out_buoy, out_beff
	logical  :: out_k, out_tke, out_ent, out_plume 
	logical	 :: out_grad, out_div, out_vort, out_dp
	logical  :: out_ccn, out_in
        logical  :: out_dtnet, out_frad, out_tmp
	logical  :: out_diagu, out_diagv, out_diagw, out_diagp
	logical  :: out_diagt, out_diagtv, out_diagq, out_diagl
	logical  :: out_diagr, out_diagi, out_diaga, out_diags, out_diagk
	logical  :: out_sat, out_flut, out_fsgs, out_prec
        logical  :: out_micro, out_aero, out_z

	logical  :: out_lwp, out_cwp, out_wvp, out_ctop, out_cmse, out_cmfl, out_mf
	logical  :: out_sfl, out_sst, out_rrate, out_srad
	logical  :: out_cape, out_cp
	logical  :: out_ints, out_olr, out_osr, out_ocs
        logical  :: out_opthic, out_zinv

	logical  :: pro_tot, pro_env, pro_cl, pro_cor
        logical  :: pro_cfev, pro_ent, pro_det
        logical  :: pro_int, pro_out, pro_halo
	logical  :: pro_up, pro_dn

        character(len=4) :: prname(10)
	integer  :: nvarout, nvarout2d, nvarsurf, nvarpro
	integer  :: n_io
	integer  :: n_av
	integer  :: n_ts
	integer  :: n_pr
	integer  :: n_sl
	integer  :: n_lag
        integer  :: n_res
        integer  :: n_pe
        integer  :: n_ne
        integer  :: n_rad
!
! --- spmd parameters:
!
      integer mypid            ! processor id
      integer nproc            ! number of processors
      integer left_p           ! left processor
      integer right_p          ! right processor
      integer front_p          ! front processor (2D decomp)
      integer back_p           ! Back processor (2D decomp)
!
! --- others
!
      logical :: new_run	 ! New run flag
      logical :: nest_run        ! Restart nested run
      logical :: rep		 ! Repeat domain at restart
      logical :: ldiag           ! Store diagnostics
      logical :: with_mom	 ! Solve for momentum
      logical :: with_scal	 ! Solve for scalars
      logical :: with_adv	 ! Include advection
      logical :: with_dif	 ! Include diffusion
      logical :: with_mic	 ! Include microphysics
      logical :: with_buoy	 ! Include buoyancy
      logical :: with_nudg       ! Include nudging
      logical :: with_lsadv      ! Include large-scale advection
      logical :: with_lssrc      ! Include large-scale sources
      logical :: with_lssub      ! Include large-scale subsidence
      logical :: with_cor        ! Include Coriolis
      logical :: with_rad        ! Include radiation
      logical :: with_tvar       ! Time dependent parameters
      logical :: with_piggy      ! Piggy-backing 
      logical :: micro_dif	 ! Turbulent diffusion for microphysics
      logical :: limit		 ! Limiter for scalar advection
      logical :: imp_buoy        ! Implicit buoyancy
      logical :: anis_k 	 ! Anisotropic turbulent diffusion 
      logical :: nl_sgs 	 ! Non-linear turbulent diffusion 
      logical :: lag_mix	 ! Particle properties from mixing with environment
      logical :: lavisc          ! Artificial viscosity
      logical :: p_mcons	 ! Mass conservation constraint on dynaical pressure
      logical :: res_lag	 ! Reset lagrangian particles when removed
      logical :: no_out          ! Turn off outputs (except last)
      logical :: minmax		 ! Write minmax values of output variables
      logical :: out_surf        ! Output surface (2D) variables
      logical :: out_hov         ! Output Hovmoller plot
      logical :: out_yav         ! Output slices averaged along y 
      logical :: lddamp        	 ! 2D Divergence damping
      logical :: conservative    ! Conservative advection
      logical :: all_rest        ! Saves all restart files
      logical :: spec_diag       ! Special diagnostics
      logical :: ldrizz		 ! Enables drizzle/rain formation 
      logical :: lvent		 ! Phase changes include ventilation
      logical :: lkin            ! Allows kinetic growth of activated droplets 
      logical :: lrime           ! Turns on riming/graupel formation 
      logical :: split_mic       ! Split microphysics
      logical :: aerosol_lag     ! Turn on lagrangial aerosol tracking with mass and condensation/evaporation
      logical :: reg_mode	 ! Regeneration mode in aerosol module     
      logical :: aero_sfc_source ! Turn on/off aerosol surface source
      logical :: ldtmic          ! Time step by microphysics scale     
      logical :: forcing_surface ! Switch of surface forcing 
      logical :: with_dropsed    ! Enables droplet sedimentation
      logical :: with_collision  ! Enables droplet collision
      logical :: with_radmod     ! Enables modification of radiation sounding (standard atmosphere)
      logical :: with_warmadv    ! Enables (warm) advection to be nudged
      logical :: with_aero_swelling ! Enables aerosol hygroscopic growth

      integer :: verbose	 ! Write debug
      integer :: ntime,nite      ! Current time step (cumulated and local)
      integer :: nstep           ! # of integration sub-steps
      integer :: ntau		 ! Number of small time-steps per step
      integer :: ldtfix          ! Fixed general time step
      integer :: isurf           ! Surface flux modelling
      integer :: momsf           ! Turns on fluxes of momentum at the surface
      integer :: lmicro          ! Microphysics level
      integer :: laero           ! Aerosol model level
      integer :: lndrop          ! Fixed or varying droplet number concentration      
      integer :: lfreeze         ! Selects freezing of drops
      integer :: looknuc	 ! Presence of lookup tables for ice nucleation
      integer :: ktrop           ! Tropopause location
      integer :: kpert		 ! Top of perturbed layer
      integer :: j_day           ! Day of the year 1-365
      integer :: sorder          ! Scheme order for integration
      integer :: mom_ord	 ! Scheme order for momentum advection
      integer :: diff_ord	 ! Scheme order for diffusion
      integer :: nsubp		 ! Number subcycles of the pressure solver
      integer :: sponge          ! Sponge layer option (1=explicit or 2=implicit)           
      integer :: rad_sw          ! Short wave radiation switch
      integer :: rad_o3          ! Interactive O3 with radiation
      integer :: lag_ord	 ! Interpolation order for lagrangian particles
      integer :: lag_init	 ! Method for initializing particles
      integer :: nlag		 ! Number of particles to be output
      integer :: nslicex	 ! Number of slices along X to output
      integer :: nslicey	 ! Number of slices along Y to output
      integer :: nslicez	 ! Number of slices along Z to output
      integer :: kout            ! Max level for output (writes levels below kout only)
      integer :: cst_cp		 ! Use constant heat capacities
      integer :: sca_set	 ! Predefined sets of scalars
      integer :: kbl		 ! Boundary layer top level
      integer :: auto		 ! Autoconversion for SB scheme
      integer :: moments	 ! Number of moments in SB scheme
      integer :: fft_solv	 ! FFT solver for 3D case
 
      character(len=5)   :: scal_adv
      character(len=3)	 :: ice_habit
      character(len=100) :: file_init
      character(len=100) :: file_output
      character(len=100) :: file_rest
      character(len=3)   :: bcl(6)
      character(len=10)  :: compos_lag       ! Name for chemical composition of particles

      real    :: ztop		 ! Top of the domain
      real    :: zbl		 ! Boundary layer height
      real    :: ztrop		 ! troposphere altitude
      real    :: z1, z2, sratio  ! Parameters of the refined area in vertical
      real    :: zl1, zl2
      real    :: tstart, tstop, time     ! Stopping and current time
      real    :: dt0             ! basic time step
      real    :: dtau            ! small time step (time-split pressure)
      real    :: dt0_i           ! Initial time step
      real    :: dt              ! Actual time step (2*dt if Leapfrog)
      real    :: limit_ts	 ! Factor limiting minimum timestep allowed (=dt0_i/limit_ts)
      real    :: iax, iav        ! output frequency for snapshots and averaging frequency
      real    :: its             ! output frequency for time-series
      real    :: ipro            ! output frequency for profiles      
      real    :: isli            ! output frequency for slices
      real    :: ires            ! output frequency for restart
      real    :: inest           ! output frequency for nesting
      real    :: ilag		 ! output frequency for parcel outputs
      real    :: iradx           ! time step of radiation in step
      real    :: cfl_max,cfl_min ! Max. allowed CFL
      real    :: lim_tol         ! Tolerance for flux limiter
      real    :: dx,dy,dz        ! intervals along x, y, & z-axis 
      real    :: xc_nest,yc_nest           ! Center of mesh
      real    :: lx_nest,ly_nest           ! Center of mesh
      real    :: psurf           ! Surface pressure
      real    :: usurf, vsurf           ! Surface velocities
      real    :: u0shift, v0shift! Translation velocities
      real    :: dpt, dqv, dw    ! init perturbation
      real    :: shf, lhf, shf0, lhf0 ! Sensible and latent heat fluxes
      real    :: ust, scf        ! Scalar and momentum flux    
      real    :: Ddiv            ! Fixed large scale divergence
      real    :: w_up            ! Fixed upwelling velocity
      real    :: sst             ! Skin surface temperature (required mainly if isurf=0)
      real    :: ssm             ! Surface relative humidity (saturation = 1)
      real    :: c_dm,c_ds       ! Surface drag coefficients
      real    :: min_w           ! Minimum background surface wind
      real    :: alb             ! Surface broadband albedo
      real    :: emi             ! Surface broadband emissivity
      real    :: zrough          ! Surface roughness
      real    :: zdec		 ! Decay height for SGS viscosity
      real    :: diff            ! Fixed diffusion coefficient
      real    :: pran            ! Turbulent Prandtl number (if <0, no sgs model for scalars)
      real    :: ice_delay       ! Microphysics delay time (in sec)
      real    :: xauto		 ! Autoconversion threshold (Seifert-Beheng)
      real    :: qauto		 ! Autoconversion threshold (Kessler)
      real    :: qthres,wthres   ! Liquid water and velocity thresholds for diagnosing clouds
      real    :: zdamp           ! Altitude above which sponge layer
      real    :: dxdamp          ! Sponge layer thickness in X
      real    :: dydamp          ! Sponge layer thickness in Y
      real    :: tdamp           ! Relaxation time scale in the sponge layer
      real    :: tnudg           ! Nudging time scale 
      real    :: tpert           ! Time of perturbations
      real    :: dtcon           ! Condensation time-step
      real    :: xnc0_d          ! Modal aerosol diameter (simple)
      real    :: xnc0_s          ! Standard deviation aerosol (simple)
      real    :: xnc0_k		 ! Mean aerosol kappa value (simple)
      real    :: xn_ccn0         ! sfc conc. of ccn or accum. sulfate
      real    :: xn_in0          ! sfc conc. of in      
      real    :: t_local         ! local time at the start of simulation
      real    :: ctr_lat         ! lat of center +-90.0
      real    :: fcor            ! Coriolis parameters (calculated)
      real    :: slicex_io(16)   ! Slice in X positions
      real    :: slicey_io(16)   ! Slice in Y positions
      real    :: slicez_io(16)   ! Slice in Z positions
      real    :: xbubble, zbubble, xrbubble		! Position initial bubble
      real    :: zrbubble, ybubble, yrbubble		! size initial bubble
      real    :: dtbubble	 ! Initial bubble perturbation
      real    :: a_lm=1.0, b_lm=0.0, c_lm=0.0, d_lm=0.0, e_lm=0.0	 ! Coefficients for linear multistep methods
      real    :: lsqt(8), lsqt0(8)
      real    :: lstem(8), lstem0(8)
      real    :: lsp(8)
      real    :: cint             ! Fraction of microphysics tendencies (for splitting)
      real    :: cdiag            ! Fraction of tendencies in diag
      real    :: mu_lag           ! First parameter of lagrangian particle distribution
      real    :: sigma_lag        ! Second parameter of lagrangian particle distribution
      real    :: des_tot, dqt_tot, drh_tot
      real    :: olr, olr_cs, ssr, toa, toa_cs, osr    ! averaged radiation diagnostics
      real    :: dtmic
      real    :: zcp = 425.
      real    :: alpha_drop       ! Shape parameter in the droplet (generalised) gamma distribution
      real    :: nu_drop          ! Scale parameter in the droplet (generalised) gamma distribution
      real    :: t_radmod         !Time point of modifying radiation sounding
      real    :: tstart1_warmadv  ! Time point in period 1 when the (warm) advection is nudged
      real    :: tstop1_warmadv   ! Time point in period 1 when the (warm) advection is terminated
      real    :: tstart2_warmadv  ! Time point in period 2 when the (warm) advection is nudged
      real    :: tstop2_warmadv   ! Time point in period 2 when the (warm) advection is terminated
      real    :: tstart3_warmadv  ! Time point in period 3 when the (warm) advection is nudged
      real    :: tstop3_warmadv   ! Time point in period 3 when the (warm) advection is terminated
      real    :: hgf             ! Hygroscopic growth factor for aerosol mixture
!
! --- input/output files
!
      character(len = 20)  :: casename
      character(len = 20)  :: gridfile
      character(len = 100) :: opdir, outpdir, prtname
      character(len = 6)   :: how_access
      character(len = 8)   :: ts_out
      character(len = 40)  :: forcing_file !Name of the forcing file to force MIMICA (sst)
      character(len = 40)  :: forcing_file_ssm !Name of the ssm forcing file
      character(len = 40)  :: forcing_file_lhf !Name of the lhf forcing file
      character(len = 40)  :: forcing_file_shf !Name of the shf forcing file
      character(len = 40)  :: warm_adv_H1 !Name of the warm advection forcing file
      character(len = 40)  :: warm_adv_dT1 !Name of the warm advection forcing file
      character(len = 40)  :: q_adv_H1 !Name of the q advection forcing file
      character(len = 40)  :: q_adv_dQ1 !Name of the q advection forcing file
      character(len = 40)  :: warm_adv_H2 !Name of the warm advection forcing file
      character(len = 40)  :: warm_adv_dT2 !Name of the warm advection forcing file
      character(len = 40)  :: q_adv_H2 !Name of the q advection forcing file
      character(len = 40)  :: q_adv_dQ2 !Name of the q advection forcing file
      character(len = 40)  :: warm_adv_H3 !Name of the warm advection forcing file
      character(len = 40)  :: warm_adv_dT3 !Name of the warm advection forcing file
      character(len = 40)  :: q_adv_H3 !Name of the q advection forcing file
      character(len = 40)  :: q_adv_dQ3 !Name of the q advection forcing file
                   
      !
      ! --- LSODE quantities
      !
      integer		   :: lrw, liw, mf
      integer, allocatable :: iwork(:)
      real, allocatable    :: rwork(:)
      real   :: facs, faci         ! Correction factors in LSODE
      
      !
      ! --- For Poisson/Helmholtz solver
      !
      integer, dimension(4) :: mgopt
      integer, dimension(:), allocatable :: iparm
      real, dimension(:), allocatable :: fparm
      real, dimension(:), allocatable :: work
      
      real(C_DOUBLE), dimension(:,:), allocatable :: wv

      integer  :: nzpres, nxpres, nypres  ! Dim. in x, y, z (ghost cells for periodicity are excluded)
      integer  :: ngx, ngy, ngz
      integer  :: iex, jey, kez
      integer  :: flagx, flagy            ! Flags for initializing tmp arrays in FFT
      integer  :: xbc, ybc, zbc           ! Flags for BC in poisson solver

      real     :: maxdiv, maxcfl, maxcflc
      real     :: dzpres
      real     :: x1pres, x2pres
      real     :: y1pres, y2pres
      real     :: z1pres, z2pres
      real     :: cx, cy
      real, dimension(nz)             :: czl, czc, czr
      real, dimension(:), allocatable :: tmpx, tmpy
            
      !
      ! --- Rate Constants:
      !
      real, dimension(nract)   :: rk
      real, dimension(nequil)  :: rke_c
      real, dimension(nequil)  :: rke_r
      real, dimension(naqreact):: rka_c
      real, dimension(naqreact):: rka_r
      !real, dimension(201:800) :: specdata ! not used
      
!
! --- one-dimensional variables
!
      real, dimension(nz) :: u00, v00, w0    ! Reference wind soundings and subsidence
      real, dimension(nz) :: qv0, qt0, rh0   ! hor. mean of qv & qt
      real, dimension(nz) :: ql0, qi0, nl0   ! hor. mean of liquid content
      real, dimension(nz) :: t0, pt0, mse0   ! hor. mean temp, pt, and moist static energy
      real, dimension(nz) :: ccn0, in0       ! Bulk aerosols
      real, dimension(nz) :: cv0             ! reference  heat capacity
      real, dimension(nz) :: k0              ! TKE
      real, dimension(nz) :: z0              ! Vertical levels
      real, dimension(nz,nscal) :: scal0	 ! Initial scalars
      real, dimension(nz) :: s_scal0
      real, dimension(nz) :: n20

      real, dimension(:,:), allocatable :: pxx      	! pressure corrections
      real, dimension(:), allocatable :: qmin, lmax, lmin, xnmin	! Hydro concentration limits

      real, dimension (:,:), allocatable :: u_nud, v_nud !,w_nud  !Wind from forcing file (SE_ATLANTIC case)
      real, dimension (:,:), allocatable :: Na_force, qt_nud      !Aerosol number concentration and total water mixing ratio 
      real, dimension (:,:), allocatable :: ta_nud, pa_force      !Temperature (SE_ATLANTIC case)
      real, dimension (:), allocatable :: ts_force, ps_force, ssm_force    !SST, SSM, and surface pressure
      real, dimension (:), allocatable :: shf_force, lhf_force  !Sensible and Latent heat flux forcing
      real, dimension (:), allocatable :: H_warm_adv1, dT_warm_adv1  !Height and dT for warm advection forcing
      real, dimension (:), allocatable :: H_q_adv1, dQ_q_adv1  !Height and dQ for moisture advection forcing
      real, dimension (:), allocatable :: H_warm_adv2, dT_warm_adv2  !Height and dT for warm advection forcing
      real, dimension (:), allocatable :: H_q_adv2, dQ_q_adv2  !Height and dQ for moisture advection forcing
      real, dimension (:), allocatable :: H_warm_adv3, dT_warm_adv3  !Height and dT for warm advection forcing
      real, dimension (:), allocatable :: H_q_adv3, dQ_q_adv3  !Height and dQ for moisture advection forcing
    
      integer :: lev_len, time_len !number of vertical levels in the forcing file

      real, dimension (nz) :: u_force,  v_force       !Wind from forcing file (SE_ATLANTIC case)
      real, dimension (nz) :: qt_force, ta_force      !Total water mixing ratio and temperature (SE_ATLANTIC case)
      real, dimension (nz) :: coef_nudge1,coef_nudge2 !Nudging coefficients (SE_ATLANTIC case)

! --- Non-uniform vertical cell sizes for advection
!      
      real 		  :: ptop	     ! Hydrostatic pressure above domain
      real                :: zmax            ! Max altitude in the domain (at nz)
      real, dimension(nz) :: fdz0            ! Cell size in z (fdz0 = dz/dz0)
      real, dimension(nz) :: dz1, dz1w, dzw       ! Centered vertical intervals
!
! --- Variables useful to store for aerosol module
!
      integer  :: bc, bcs
      integer  :: anu, ait, acc
!
! --- MPI specific variables
!
      integer :: nxt, nyt, nyzt, nxtp, nytp, nyztp
      integer :: BLOCK3D, BLOCK3Dp, BLOCK2D_r, BLOCK2D_i, BLOCK2Dp
      integer :: GHOSTX_send_s, GHOSTX_send_e, GHOSTY_send_s, GHOSTY_send_e
      integer :: GHOSTX_send_s_2d, GHOSTX_send_e_2d, GHOSTY_send_s_2d, GHOSTY_send_e_2d
        
      end module shared_data

