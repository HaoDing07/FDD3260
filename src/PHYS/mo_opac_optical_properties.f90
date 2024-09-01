MODULE mo_salsa_optical_properties
	!USE classSpecies, ONLY : maxspec
	!USE mo_submctl, ONLY : spec
	IMPLICIT NONE
      
! Refractive indices and the corresponding wavelengths
! Shortwave
INTEGER, PARAMETER :: NaerRadPropsSW = 23 
REAL, PARAMETER :: aerRefrIBands_SW(NaerRadPropsSW) =  (/ 0.25,  0.3 ,  0.35,  0.4 ,  0.45,  0.5 ,  0.55,  0.6 ,  0.65,       &
								0.7 ,  0.75,  0.8 ,  0.9 ,  1.  ,  1.25,  1.5 ,  1.75,  2.  , &
								2.5 ,  3.  ,  3.2 ,  3.39,  3.5 /) *1.e-4  ! cm

REAL, PARAMETER :: riReSWsu(NaerRadPropsSW) =                                                                             &
						       (/1.484, 1.469, 1.452, 1.44 , 1.432, 1.431, 1.43 , 1.429, 1.429,   &
							1.428, 1.427, 1.426, 1.425, 1.422, 1.413, 1.403, 1.394, 1.384,    &
							1.344, 1.293, 1.311, 1.35 , 1.376/),                              &

		   riImSWsu(NaerRadPropsSW) =                                                                             & 
		   				      (/1.000e-08, 1.000e-08, 1.000e-08, 1.000e-08, 1.000e-08, 1.000e-08, &
		   					1.000e-08, 1.470e-08, 1.670e-08, 2.050e-08, 7.170e-08, 8.630e-08, &
		   					2.550e-07, 1.530e-06, 6.940e-06, 1.200e-04, 4.160e-04, 1.260e-03, &
		   					3.760e-03, 9.550e-02, 1.350e-01, 1.578e-01, 1.580e-01/)    							   
		
REAL, PARAMETER :: riReSWbc(NaerRadPropsSW) =                                                                              &
						      (/1.62, 1.74, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75,  &
							1.75, 1.75, 1.76, 1.76, 1.77, 1.79, 1.8 , 1.82, 1.84, 1.86, 1.87,  &
							1.88/),                                                            &
		   riImSWbc(NaerRadPropsSW) =                                                                               &
						      (/0.45  , 0.47  , 0.465 , 0.46  , 0.455 , 0.45  , 0.44  , 0.435 ,     &
							0.435 , 0.43  , 0.43  , 0.43  , 0.435 , 0.44  , 0.45  , 0.46  ,     &
							0.48  , 0.49  , 0.51  , 0.54  , 0.54  , 0.5495, 0.56  /)            

REAL, PARAMETER :: riReSWoc(NaerRadPropsSW) =									     &
						(/1.53 , 1.53 , 1.53 , 1.53 , 1.53 , 1.53 , 1.53 , 1.53 , 1.53 ,     &
						1.53 , 1.53 , 1.52 , 1.52 , 1.52 , 1.46 , 1.41 , 1.34 , 1.26 ,       &
						1.18 , 1.16 , 1.22 , 1.258, 1.28 /),                                 &
		   riImSWoc(NaerRadPropsSW) =                                                                        &
						(/0.03   , 0.008  , 0.008  , 0.008  , 0.008  , 0.008  , 0.008  ,     &
						0.008  , 0.008  , 0.008  , 0.008  , 0.008  , 0.008  , 0.008  ,       &
						0.008  , 0.008  , 0.008  , 0.008  , 0.009  , 0.012  , 0.01   ,       &
						0.01285, 0.011  /)                                                   

REAL, PARAMETER :: riReSWss(NaerRadPropsSW) =                                                                         &
						(/1.51, 1.51, 1.51, 1.5 , 1.5 , 1.5 , 1.5 , 1.49, 1.49, 1.49, 1.49,   &
						  1.48, 1.48, 1.47, 1.47, 1.46, 1.45, 1.45, 1.43, 1.61, 1.49, 1.48,   &
						  1.48/),                                                             &						
		   riImSWss(NaerRadPropsSW) =                                                                          &
						  (/5.00e-06, 2.00e-06, 3.24e-07, 3.00e-08, 2.43e-08, 1.55e-08,        &
						    1.00e-08, 1.60e-08, 4.24e-08, 2.00e-07, 1.08e-06, 1.95e-06,        &
						    4.24e-05, 1.41e-04, 3.58e-04, 5.70e-04, 7.62e-04, 1.00e-03,        &
						    4.00e-03, 1.00e-02, 3.00e-03, 2.05e-03, 1.60e-03/)
                 								
						
! Longwave							
INTEGER, PARAMETER :: NaerRadPropsLW = 38
REAL, PARAMETER :: aerRefrIBands_LW(NaerRadPropsLW) =                                                                  & 
							(/3.75,  4.  ,  4.5 ,  5.  ,                                   &
							5.5 ,  6.  ,  6.2 ,  6.5 ,  7.2 ,  7.9 ,  8.2 ,  8.5 ,  8.7 ,  &
							9.  ,  9.2 ,  9.5 ,  9.8 , 10.  , 10.6 , 11.  , 11.5 , 12.5 ,  &
							13.  , 14.  , 14.8 , 15.  , 16.4 , 17.2 , 18.  , 18.5 , 20.  , &
							21.3 , 22.5 , 25.  , 27.9 , 30.  , 35.  , 40.  /) *1.e-4  ! cm

REAL, PARAMETER :: riReLWsu(NaerRadPropsLW) =                                                                             &
							(/1.396, 1.385, 1.385, 1.36 , 1.337, 1.425, 1.424, 1.37 , 1.21 ,  &
							1.14 , 1.2  , 1.37 , 1.53 , 1.65 , 1.6  , 1.67 , 1.91 , 1.89 ,    &
							1.72 , 1.67 , 1.89 , 1.74 , 1.69 , 1.64 , 1.61 , 1.59 , 1.52 ,    &
							1.724, 1.95 , 1.927, 1.823, 1.78 , 1.87 , 1.93 , 1.92 , 1.92 ,    &
							1.9  , 1.89 /),                                                   &
		   riImLWsu(NaerRadPropsLW) =                                                                              &
							(/0.131  , 0.126  , 0.12   , 0.121  , 0.183  , 0.195  , 0.165  ,   &
							0.128  , 0.176  , 0.488  , 0.645  , 0.755  , 0.772  , 0.633  ,     &
							0.586  , 0.75   , 0.68   , 0.455  , 0.34   , 0.485  , 0.374  ,     &
							0.198  , 0.195  , 0.195  , 0.205  , 0.211  , 0.414  , 0.59   ,     &
							0.041  , 0.03025, 0.02352, 0.2925 , 0.0315 , 0.2    , 0.18   ,     &
							0.18   , 0.19   , 0.22   /)                                        

REAL, PARAMETER :: riReLWbc(NaerRadPropsLW) =                                                                               &
							(/1.9 , 1.92, 1.94, 1.97, 1.99, 2.02, 2.03, 2.04, 2.06, 2.12, 2.13, &
							2.15, 2.16, 2.17, 2.18, 2.19, 2.2 , 2.21, 2.22, 2.23, 2.24, 2.27,   &
							2.28, 2.31, 2.33, 2.33, 2.36, 2.38, 2.4 , 2.41, 2.45, 2.46, 2.48,   &
							2.51, 2.54, 2.57, 2.63, 2.69/),                                     &
		   riImLWbc(NaerRadPropsLW) =                                                                             &
							(/0.57 , 0.58 , 0.59 , 0.6  , 0.61 , 0.62 , 0.625, 0.63 , 0.65 ,  &
							0.67 , 0.68 , 0.69 , 0.69 , 0.7  , 0.7  , 0.71 , 0.715, 0.72 ,    &
							0.73 , 0.73 , 0.74 , 0.75 , 0.76 , 0.775, 0.79 , 0.79 , 0.81 ,    &
							0.82 , 0.825, 0.83 , 0.85 , 0.86 , 0.87 , 0.89 , 0.91 , 0.93 ,    &
							0.97 , 1.   /)                                                     
							
REAL, PARAMETER :: riReLWoc(NaerRadPropsLW) =                                                                             &
						      (/1.27, 1.26, 1.26, 1.25, 1.22, 1.15, 1.14, 1.13, 1.4 , 1.15, 1.13, &
							1.3 , 1.4 , 1.7 , 1.72, 1.73, 1.74, 1.75, 1.62, 1.62, 1.59, 1.51, &
							1.47, 1.52, 1.57, 1.57, 1.6 , 1.63, 1.64, 1.64, 1.68, 1.77, 1.9 , &
							1.97, 1.89, 1.8 , 1.9 , 2.1 /),                                   &
		   riImLWoc(NaerRadPropsLW) =                                                                             &
						      (/0.011 , 0.012 , 0.014 , 0.016 , 0.021 , 0.037 , 0.039 , 0.042 ,   & 
							0.055 , 0.04  , 0.0742, 0.09  , 0.1   , 0.14  , 0.15  , 0.162 ,   &
							0.162 , 0.162 , 0.12  , 0.105 , 0.1   , 0.09  , 0.1   , 0.085 ,   &
							0.1   , 0.1   , 0.1   , 0.1   , 0.115 , 0.12  , 0.22  , 0.28  ,   &
							0.28  , 0.24  , 0.32  , 0.42  , 0.5   , 0.6   /)		  					  

REAL, PARAMETER :: riReLWss(NaerRadPropsLW) =                                                                                &
							(/1.47, 1.48, 1.49, 1.47, 1.42, 1.41, 1.6 , 1.46, 1.42, 1.4 , 1.42,  &
							1.48, 1.6 , 1.65, 1.61, 1.58, 1.56, 1.54, 1.5 , 1.48, 1.48, 1.42,    &
							1.41, 1.41, 1.43, 1.45, 1.56, 1.74, 1.78, 1.77, 1.76, 1.76, 1.76,    &
							1.76, 1.77, 1.77, 1.76, 1.74/),                                      &
		   riImLWss(NaerRadPropsLW) =                                                                                &
							(/0.0014, 0.0014, 0.0014, 0.0025, 0.0036, 0.011 , 0.022 , 0.005 ,    &
							0.007 , 0.013 , 0.02  , 0.026 , 0.03  , 0.028 , 0.0262, 0.018 ,      &
							0.016 , 0.015 , 0.014 , 0.014 , 0.014 , 0.016 , 0.018 , 0.023 ,      &
							0.03  , 0.035 , 0.09  , 0.12  , 0.13  , 0.135 , 0.152 , 0.165 ,      &
							0.18  , 0.205 , 0.275 , 0.3   , 0.5   , 1.    /)    
                 
      
	! allocatable arrays for the chemicals configured active in SALSA (order as in the NAMELIST and mass arrays)
	! Shape = (number of species, number of spectral bands)
	REAL, ALLOCATABLE :: riReSW(:,:), riImSW(:,:), riReLW(:,:), riImLW(:,:)
      
	CONTAINS
	    
	SUBROUTINE initialize_optical_properties()
	  IMPLICIT NONE
	  
	  INTEGER :: nspec, ii, jj
      
	  !###################  This part does not appear in the original file in UCLALES-SALSA
      
	  !This is taken from src/src_shared/classSpecies.f90 ========================
	 INTEGER,PARAMETER :: maxspec = 4  ! Maximum number of aerosol species, excluding water    
      
	  !This appears in mo_salsa_init.f90 =================================:
	  CHARACTER(len=4), TARGET   :: allNames(maxspec)  = ['SO4 ','BC  ','SS','OC']   
      
	  INTEGER                    :: allInd(maxspec) = [1,2,3,4]       ! Indices of compounds in model mass arrays (0 if not active)
      
      
      
      
	  ! Refractive indices for all available chemicals (order as in classSpecies/allNames) 
	  REAL :: riReSW_all(maxspec,NaerRadPropsSW), riImSW_all(maxspec,NaerRadPropsSW),   &
		  riReLW_all(maxspec,NaerRadPropsLW), riImLW_all(maxspec,NaerRadPropsLW)
	  
	  ! The number of active compounds defined in the NAMELIST (+water)
	  !nspec = spec%getNSpec(type="wet")
      
	  !Testing a number for nspec
	  nspec = 4                           
	  !======================
      
	  ! This gives the optical properties for all possible SALSA compounds in single arrays, hard coded in the same order
	  ! as classSpecies/allNames. This is mostly needed because it makes it easy to use the indexing provided by the class object "spec" to
	  ! construct the trucated arrays below
	  riReSW_all(:,:) = RESHAPE( [riReSWsu,riReSWbc,riReSWss,riReSWoc], &
				     SHAPE(riReSW_all), ORDER=[2,1]   )
	  riReLW_all(:,:) = RESHAPE( [riReLWsu,riReLWbc,riReLWss,riReLWoc], &
				     SHAPE(riReLW_all), ORDER=[2,1]   )
	  riImSW_all(:,:) = RESHAPE( [riImSWsu,riImSWbc,riImSWss,riImSWoc], &
				     SHAPE(riImSW_all), ORDER=[2,1]   )
	  riImLW_all(:,:) = RESHAPE( [riImLWsu,riImLWbc,riImLWss,riImLWoc], &
				     SHAPE(riImLW_all), ORDER=[2,1]   )
			       
	  ! Allocate and construct the truncated optical property arrays (contain only active compounds and given in the order in which they
	  ! are given in the NAMELIST and the LES/SALSA mass arrays)
	  ALLOCATE( riReSW(nspec,NaerRadPropsSW), riImSW(nspec,NaerRadPropsSW),  &
		    riReLW(nspec,NaerRadPropsLW), riImLW(nspec,NaerRadPropsLW)   )
	  
	  riReSW = 0.; riImSW = 0.; riReLW = 0.; riImLW = 0.
      
	  ! Sort the properties
	  DO ii = 1,maxspec
	     !jj = spec%allInd(ii)
	      jj =allInd(ii)
	     IF (jj == 0) CYCLE  ! Species not used
	     ! DEBUGGING
	     IF (jj > nspec) STOP "optical props vaarin meni!!!"
	     riReSW(jj,:) = riReSW_all(ii,:)
	     riImSW(jj,:) = riImSW_all(ii,:)
	     riReLW(jj,:) = riReLW_all(ii,:)
	     riImLW(jj,:) = riImLW_all(ii,:)
	  END DO
	     
	  print *,"Optical properties initialized"
      
	END SUBROUTINE initialize_optical_properties
      
      
      END MODULE mo_salsa_optical_properties


