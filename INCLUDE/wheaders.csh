#!/bin/csh -f

rm $INCDIR/ctrparam.h
rm ./compile.log

cat >>! ./compile.log << EOF
Compiling MIMICA V5 version $VERSION
EOF
cat >>! $INCDIR/ctrparam.h << EOF
#define VERSION $VERSION
EOF

if ($?CMPLER) then
if ($CMPLER == "INTEL") then
cat >>! ./compile.log << EOF
Compiling with the intel compiler
EOF
endif
endif

if ($DEBUG == "TRUE") then
cat >>! ./compile.log << EOF
Compiling with debug flags
EOF
else if ($PROF == "TRUE") then
cat >>! ./compile.log << EOF
Compiling with profiling flags
EOF
else
cat >>! ./compile.log << EOF
Compiling with optimization flags
EOF
endif

cat >>! ./compile.log << EOF

Summary of options selected in "start":
EOF

if ($?DEBUG) then
if ($DEBUG == "TRUE") then
cat >>! ./compile.log << EOF
#define DEBUG
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef DEBUG
EOF
cat >>! ./compile.log << EOF
#undefine DEBUG
EOF
endif
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef DEBUG
EOF
cat >>! ./compile.log << EOF
Missing option keyword: DEBUG 
EOF
endif

if ($?PROF) then
if ($PROF == "TRUE") then
cat >>! ./compile.log << EOF
#define PROF
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef PROF
EOF
cat >>! ./compile.log << EOF
#undefine PROF
EOF
endif
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef PROF
EOF
cat >>! ./compile.log << EOF
Missing option keyword: PROF 
EOF
endif

if ($?SPMD) then
if ($SPMD == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define SPMD
#define INTTYPE  MPI_INTEGER
#define LOGTYPE  MPI_LOGICAL
#define CARTYPE  MPI_CHARACTER
#define REALTYPE MPI_REAL8
EOF
cat >>! ./compile.log << EOF
#define SPMD
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef SPMD
EOF
cat >>! ./compile.log << EOF
#undefine SPMD
EOF
endif
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef SPMD
EOF
cat >>! ./compile.log << EOF
Missing option keyword: PROF 
EOF
endif

# ===
# === choose either 3d or 2d  and domain decomposition
# ===

if ($?M3D) then
if ($M3D == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define MODEL_3D 1
EOF
cat >>! ./compile.log << EOF
#define M3D
EOF
else 
cat >>! $INCDIR/ctrparam.h << EOF
#undef MODEL_3D
EOF
cat >>! ./compile.log << EOF
#undefine M3D
EOF
endif
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef MODEL_3D
EOF
cat >>! ./compile.log << EOF
Missing option keyword: M3D 
EOF
endif

if ($?COLUMN) then
if ($COLUMN == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define COLUMN 1
EOF
cat >>! ./compile.log << EOF
#define COLUMN
EOF
else 
cat >>! $INCDIR/ctrparam.h << EOF
#undef COLUMN
EOF
cat >>! ./compile.log << EOF
#undefine COLUMN
EOF
endif
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef COLUMN
EOF
cat >>! ./compile.log << EOF
Missing option keyword: COLUMN 
EOF
endif

if ($?DECOMP_2D) then
if ($DECOMP_2D == "TRUE" && $SPMD == "TRUE" && $M3D == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define DECOMP_2D 1
EOF
cat >>! ./compile.log << EOF
#define DECOMP_2D
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef DECOMP_2D
EOF
cat >>! ./compile.log << EOF
#undefine DECOMP_2D
EOF
endif
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef DECOMP_2D
EOF
cat >>! ./compile.log << EOF
Missing option keyword: DECOMP_2D 
EOF
endif

cat >>! ./compile.log << EOF
#MAXX $MAXX
#MAXY $MAXY
#MAXZ $MAXZ
#NPX $NPX
#NPY $NPY
EOF

if ($?FINE) then
if ($FINE == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define FINE 1
EOF
cat >>! ./compile.log << EOF
#define FINE
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef FINE
EOF
cat >>! ./compile.log << EOF
#undefine FINE
EOF
endif
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef FINE
EOF
cat >>! ./compile.log << EOF
Missing option keyword: FINE 
EOF
endif

if ($?NESTING) then
if ($NESTING == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define NESTING 1
EOF
cat >>! ./compile.log << EOF
#define NESTING
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef NESTING
EOF
cat >>! ./compile.log << EOF
#undefine NESTING
EOF
endif
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef NESTING
EOF
cat >>! ./compile.log << EOF
Missing option keyword: NESTING 
EOF
endif

if ($?CHANNEL) then
if ($CHANNEL == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define CHANNEL 1
EOF
cat >>! ./compile.log << EOF
#define CHANNEL
EOF
else 
cat >>! $INCDIR/ctrparam.h << EOF
#undef CHANNEL
EOF
cat >>! ./compile.log << EOF
#undefine CHANNEL
EOF
endif
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef CHANNEL
EOF
cat >>! ./compile.log << EOF
Missing option keyword: CHANNEL 
EOF
endif

if ($?PARALLEL_OUT) then
if ($PARALLEL_OUT == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define PARALLEL_OUT 1
EOF
cat >>! ./compile.log << EOF
#define PARALLEL_OUT
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef PARALLEL_OUT
EOF
cat >>! ./compile.log << EOF
#undefine PARALLEL_OUT
EOF
endif
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef PARALLEL_OUT
EOF
cat >>! ./compile.log << EOF
Missing option keyword: PARALLEL_OUT 
EOF
endif

if ($?MEANDATA) then
if ($MEANDATA == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define MEANDATA 1
EOF
cat >>! ./compile.log << EOF
#define MEANDATA
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef MEANDATA
EOF
cat >>! ./compile.log << EOF
#undefine MEANDATA
EOF
endif
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef MEANDATA
EOF
cat >>! ./compile.log << EOF
Missing option keyword: MEANDATA 
EOF
endif

if ($?CLOUD_TRACK) then
cat >>! ./compile.log << EOF
Deprecated option keyword: CLOUD_TRACK 
EOF
endif

if ($?CORIOLIS) then
cat >>! ./compile.log << EOF
Deprecated option keyword: CORIOLIS 
EOF
endif

if ($?RADIA) then
if ($RADIA == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define RAD_ENABLE 1
EOF
cat >>! ./compile.log << EOF
#define RADIA
EOF
else 
cat >>! $INCDIR/ctrparam.h << EOF
#undef RAD_ENABLE
EOF
cat >>! ./compile.log << EOF
#undefine RADIA
EOF
endif
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef RAD_ENABLE
EOF
cat >>! ./compile.log << EOF
Missing option keyword: RADIA 
EOF
endif

if ($?CHEM) then
if ($CHEM == "TRUE") then
cat >>! ./compile.log << EOF
#define CHEM
EOF
else
cat >>! ./compile.log << EOF
#undefine CHEM
EOF
endif

if ($?AEROSOL) then
if ($AEROSOL == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define AERO_ENABLE
EOF
cat >>! ./compile.log << EOF
#define AEROSOL
EOF
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef AERO_ENABLE
EOF
cat >>! ./compile.log << EOF
#undefine AEROSOL
EOF
endif
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef AEROSOL
EOF
cat >>! ./compile.log << EOF
Missing option keyword: AEROSOL 
EOF
endif
endif

if ($?ADV_AEROCHEM) then
cat >>! ./compile.log << EOF
Deprecated option keyword: ADV_AEROCHEM
EOF
endif

if ($?LSADV) then
cat >>! ./compile.log << EOF
Deprecated option keyword: LSADV 
EOF
endif

if ($?TKE) then
if ($TKE == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define TKE 1
EOF
cat >>! ./compile.log << EOF
#define TKE
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef TKE
EOF
cat >>! ./compile.log << EOF
#undefine TKE
EOF
endif
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef TKE
EOF
cat >>! ./compile.log << EOF
Missing option keyword: TKE 
EOF
endif

if ($?LAGRANGE) then
if ($LAGRANGE == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define LAGRANGE 1
EOF
cat >>! ./compile.log << EOF
#define LAGRANGE
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef LAGRANGE
EOF
cat >>! ./compile.log << EOF
#undefine LAGRANGE
EOF
endif
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef LAGRANGE
EOF
cat >>! ./compile.log << EOF
Missing option keyword: LAGRANGE 
EOF
endif

if ($?ISENTROPIC) then
if ($ISENTROPIC == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define ISENTROPIC 1
EOF
cat >>! ./compile.log << EOF
#define ISENTROPIC
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef ISENTROPIC
EOF
cat >>! ./compile.log << EOF
#undefine ISENTROPIC
EOF
endif
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef ISENTROPIC
EOF
cat >>! ./compile.log << EOF
Missing option keyword: ISENTROPIC 
EOF
endif

if ($?ANELASTIC) then
if ($ANELASTIC == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define ANELASTIC 1
EOF
cat >>! ./compile.log << EOF
#define ANELASTIC
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef ANELASTIC
EOF
cat >>! ./compile.log << EOF
#undefine ANELASTIC
EOF
endif
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef ANELASTIC
EOF
cat >>! ./compile.log << EOF
Missing option keyword: ANELASTIC 
EOF
endif

if ($?FFT_SOLVER) then
cat >>! ./compile.log << EOF
Deprecated option keyword: FFT_SOLVER 
EOF
endif

if ($?CONSERVATIVE) then
if ($CONSERVATIVE == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define CONSERVATIVE 1
EOF
cat >>! ./compile.log << EOF
#define CONSERVATIVE
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef CONSERVATIVE
EOF
cat >>! ./compile.log << EOF
#undefine CONSERVATIVE
EOF
endif
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef CONSERVATIVE
EOF
cat >>! ./compile.log << EOF
Missing option keyword: CONSERVATIVE 
EOF
endif

if ($?ADV_SPLIT) then
if ($ADV_SPLIT == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define ADV_SPLIT 1
EOF
cat >>! ./compile.log << EOF
#define ADV_SPLIT
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef ADV_SPLIT
EOF
cat >>! ./compile.log << EOF
#undefine ADV_SPLIT
EOF
endif
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef ADV_SPLIT
EOF
cat >>! ./compile.log << EOF
Missing option keyword: ADV_SPLIT 
EOF
endif

if ($?RK) then
cat >>! ./compile.log << EOF
Deprecated option: RK
EOF
endif

if ($?SAT_ADJ) then
if ($SAT_ADJ == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define SAT_ADJ 1
EOF
cat >>! ./compile.log << EOF
#define SAT_ADJ
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef SAT_ADJ
EOF
cat >>! ./compile.log << EOF
#undefine SAT_ADJ
EOF
endif
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef SAT_ADJ
EOF
cat >>! ./compile.log << EOF
Missing option keyword: SAT_ADJ 
EOF
endif

if ($?SEIFERT) then
if ($SEIFERT == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define SEIFERT 1
EOF
cat >>! ./compile.log << EOF
#define SEIFERT
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef SEIFERT
EOF
cat >>! ./compile.log << EOF
#undefine SEIFERT
EOF
endif
else
setenv SEIFERT FALSE
cat >> ! $INCDIR/ctrparam.h << EOF
#undef SEIFERT
EOF
cat >>! ./compile.log << EOF
Missing option keyword: SEIFERT 
EOF
endif

if ($?NUC_CNT) then
if ($NUC_CNT == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define NUC_CNT 1
EOF
cat >>! ./compile.log << EOF
#define NUC_CNT
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef NUC_CNT
EOF
cat >>! ./compile.log << EOF
#undefine NUC_CNT
EOF
endif
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef NUC_CNT
EOF
cat >>! ./compile.log << EOF
Missing option keyword: NUC_CNT 
EOF
endif

cat >> ! $INCDIR/ctrparam.h << EOF
#undef NUC_CNT1
EOF
cat >>! ./compile.log << EOF
#undefine NUC_CNT1
EOF

if ($?AERO_SCAV) then
cat >>! ./compile.log << EOF
Deprecated option keyword: AERO_SCAV 
EOF
endif

if ($?AERO_RADIA) then
if ($AERO_RADIA == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF
#define AERO_RADIA 1
EOF
cat >>! ./compile.log << EOF
#define AERO_RADIA
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef AERO_RADIA
EOF
cat >>! ./compile.log << EOF
#undefine AERO_RADIA
EOF
endif
else
cat >> ! $INCDIR/ctrparam.h << EOF
#undef AERO_RADIA
EOF
cat >>! ./compile.log << EOF
Missing option keyword: AERO_RADIA 
EOF
 endif



! === 
! === Set options specific to AEROSOL/CHEM
! ===

cat >>! $INCDIR/ctrparam.h << EOF

EOF
cat >>! ./compile.log << EOF

EOF

cat >>! $INCDIR/ctrparam.h << EOF
#define NVRAD $MAXZ
EOF

if ($CHEM == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF

! ===
! === if include chemistry 
! === (either transport only or transport + reaction)
! ===
!#define CHEM_ENABLE 1

! ===
! === if include chemical reactions
! ===
#define CHEM_REACT 1
!#undef CHEM_REACT
! ===
! === if include aqueous chemical reactions
! ===
#define AQCHEM_ENABLE 1
!#undef AQCHEM_ENABLE
! ===
! === if include heterogeneous chemical reactions
! ===
#define SOLIDCHEM_ENABLE 1
!#undef SOLIDCHEM_ENABLE
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF

#undef CHEM_ENABLE
#undef CHEM_REACT
#undef AQCHEM_ENABLE
#undef SOLIDCHEM_ENABLE
EOF
endif

if ($AEROSOL == "TRUE") then
cat >>! $INCDIR/ctrparam.h << EOF

! ===
! === if calculate aerosol distribution 
! ===
#define AERO_ENABLE 
! ===
! === if calculate aerosol condensation 
! ===
!#define AERO_COND 
#undef AERO_COND
! ===
! === if calculate aerosol coagulation 
! ===
!#define AERO_COAG 
#undef AERO_COAG
! ===
! === if calculate aerosol dry deposition
! ===
!#undef AERO_DDEP
#undef AERO_DDEP 
! ===
! === if check aerosol size interval
! ===
!#define AERO_CHECK 
#undef AERO_CHECK
! ===
! === if ageing of pure carbon aerosols 
! ===
!#define AERO_AGE 
#undef AERO_AGE
! ===
! === if calculate aerosol nucleation 
! ===
#undef AERO_NUC 
! ===
! === if new binary nucleation   
! ===
!#undef AERO_NUCBIN 
! ===
! === if new ternary nucleation   
! ===
!#undef AERO_NUCTERN 
! ===
! === if ammonia    
! ===
#undef AMMONIA 
! ===
EOF
else
cat >>! $INCDIR/ctrparam.h << EOF
#undef AERO_ENABLE
#undef IMPACT_OLD
EOF
endif
