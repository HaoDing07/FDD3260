
! ================================================================
!
!  Revision:
!	Date	By		Brief Description
!	----	--		-----------------	
!	041902	Chien Wang	add heterogeneous chemistry
!
! ================================================================

! --- define chemical variables

	! --- mass mixing ratios of gaseous chemicals
#ifdef CHEM_ENABLE
#if ( defined SPMD )
	type (gas_chemical),                           &
       dimension(ip_start:ip_end,jp_start:jp_end,1:nz) & 
     				:: gas
#else
        type (gas_chemical), dimension(nx,ny,nz) :: gas
#endif
#else
	type (gas_chemical) :: gas
#endif

	! --- mass mixing ratios of aqueous chemicals
#ifdef AQCHEM_ENABLE
#if ( defined SPMD )
	type (aq_chemical),                             & 
       dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  &
     				:: aqc, aqr
#else
        type (aq_chemical),  dimension(nx,ny,nz) :: aqc, aqr
#endif
#else
	type (gas_chemical) :: aqc, aqr
#endif


	! --- mass mixing ratios of solid phase chemicals
#ifdef SOLIDCHEM_ENABLE
#if ( defined SPMD )
	type (solid_chemical),                           & 
       dimension(ip_start:ip_end,jp_start:jp_end,1:nz)   & 
     				:: solidi
#else
        type (solid_chemical),  dimension(nx,ny,nz) :: solidi
#endif
#else
	type (gas_chemical) :: solidi
#endif

