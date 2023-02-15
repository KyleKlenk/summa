! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module systemSolv_module
USE, intrinsic :: iso_c_binding
! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access matrix information
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! global metadata
USE globalData,only:flux_meta                        ! metadata on the model fluxes

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J K-1)
                    Tfreeze,      & ! temperature at freezing              (K)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookFORCE      ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:       &
 localColumn,                     & ! separate groundwater representation in each local soil column
 singleBasin                        ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:      &
 qbaseTopmodel,                  & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                      & ! a big bucket (lumped aquifer model)
 noExplicit                        ! no explicit groundwater parameterization

! safety: set private unless specified otherwise
implicit none
private
public::systemSolv
public::free_kinsol_objects

! control parameters
real(rkind),parameter  :: valueMissing=-9999._rkind     ! missing value
real(rkind),parameter  :: verySmall=1.e-12_rkind        ! a very small number (used to check consistency)
real(rkind),parameter  :: veryBig=1.e+20_rkind          ! a very big number
real(rkind),parameter  :: dx = 1.e-8_rkind              ! finite difference increment

! Sundials Parameters
real(c_double),  parameter :: ftol = 1.d-5
real(c_double),  parameter :: stol = 1.d-5

contains


! **********************************************************************************************************
! public subroutine systemSolv: run the coupled energy-mass model for one timestep
! **********************************************************************************************************
subroutine systemSolv(&
                       ! input: model control
                       dt,                & ! intent(in):    time step (s)
                       nState,            & ! intent(in):    total number of state variables
                       firstSubStep,      & ! intent(in):    flag to denote first sub-step
                       firstFluxCall,     & ! intent(inout): flag to indicate if we are processing the first flux call
                       firstSplitOper,    & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                       computeVegFlux,    & ! intent(in):    flag to denote if computing energy flux over vegetation
                       scalarSolution,    & ! intent(in):    flag to denote if implementing the scalar solution
                       ! input/output: data structures
                       type_data,         & ! intent(in):    type of vegetation and soil
                       attr_data,         & ! intent(in):    spatial attributes
                       forc_data,         & ! intent(in):    model forcing data
                       mpar_data,         & ! intent(in):    model parameters
                       indx_data,         & ! intent(inout): index data
                       prog_data,         & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,         & ! intent(inout): model diagnostic variables for a local HRU
                       flux_temp,         & ! intent(inout): model fluxes for a local HRU
                       bvar_data,         & ! intent(in):    model variables for the local basin
                       model_decisions,   & ! intent(in):    model decisions
                       stateVecInit,      & ! intent(in):    initial state vector
                       ! output
                       deriv_data,        & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       ixSaturation,      & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                       untappedMelt,      & ! intent(out):   un-tapped melt energy (J m-3 s-1)
                       stateVecTrial,     & ! intent(out):   updated state vector
                       reduceCoupledStep, & ! intent(out):   flag to reduce the length of the coupled step
                       tooMuchMelt,       & ! intent(out):   flag to denote that there was too much melt
                       niter,             & ! intent(out):   number of iterations taken
                       err,message)         ! intent(out):   error code and error message
  ! ---------------------------------------------------------------------------------------
  ! structure allocations
  USE allocspace_module,only:allocLocal                ! allocate local data structures
  ! simulation of fluxes and residuals given a trial state vector
  USE eval8summa_module,only:eval8summa                ! simulation of fluxes and residuals given a trial state vector
  USE eval8summa_module,only:eval8summa_kinsol
  USE summaSolve_module,only:summaSolve                ! calculate the iteration increment, evaluate the new state, and refine if necessary
  USE getVectorz_module,only:getScaling                ! get the scaling vectors
  USE convE2Temp_module,only:temp2ethpy                ! convert temperature to enthalpy
  USE computJacob_module,only:computJacob_kinsol
  
  ! Sundials modules
  USE fsundials_context_mod                            ! Fortran interface to SUNContext
  USE fkinsol_mod                                      ! Fortran interface to KINSOL
  USE fnvector_serial_mod                              ! Fortran interface to N_Vector
  USE fsundials_nvector_mod                            ! Fortran interface to SUNDIALS N_Vector
  USE fsunmatrix_dense_mod                             ! Fortran interface to SUNMatrix
  USE fsundials_matrix_mod                             ! Fortran interface to SUNDIALS Matrix
  USE fsunlinsol_dense_mod                             ! Fortran interface to dense SUNLinearSolver
  USE fsundials_linearsolver_mod                       ! Fortran interface to generic SUNLinearSolver
  
  USE kinsol_user_data_type                            ! user data type for KINSOL
 USE var_lookup,only:maxvarDecisions ! maximum number of decisions
  
  implicit none
  ! ---------------------------------------------------------------------------------------
  ! * dummy variables
  ! ---------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)          :: dt                            ! time step (seconds)
  integer(i4b),intent(in)         :: nState                        ! total number of state variables
  logical(lgt),intent(in)         :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(inout)      :: firstFluxCall                 ! flag to define the first flux call
  logical(lgt),intent(in)         :: firstSplitOper                ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt),intent(in)         :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  logical(lgt),intent(in)         :: scalarSolution                ! flag to denote if implementing the scalar solution
  ! input/output: data structures
  type(var_i),intent(in)          :: type_data                     ! type of vegetation and soil
  type(var_d),intent(in)          :: attr_data                     ! spatial attributes
  type(var_d),intent(in)          :: forc_data                     ! model forcing data
  type(var_dlength),intent(in)    :: mpar_data                     ! model parameters
  type(var_ilength),intent(inout) :: indx_data                     ! indices for a local HRU
  type(var_dlength),intent(inout) :: prog_data                     ! prognostic variables for a local HRU
  type(var_dlength),intent(inout) :: diag_data                     ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_temp                     ! model fluxes for a local HRU
  type(var_dlength),intent(in)    :: bvar_data                     ! model variables for the local basin
  type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
  real(rkind),intent(in)          :: stateVecInit(:)               ! initial state vector (mixed units)
  ! output: model control
  type(var_dlength),intent(inout) :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
  integer(i4b),intent(inout)      :: ixSaturation                  ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
  real(rkind),intent(out)         :: untappedMelt(:)               ! un-tapped melt energy (J m-3 s-1)
  real(rkind),intent(out)         :: stateVecTrial(:)              ! trial state vector (mixed units)
  logical(lgt),intent(out)        :: reduceCoupledStep             ! flag to reduce the length of the coupled step
  logical(lgt),intent(out)        :: tooMuchMelt                   ! flag to denote that there was too much melt
  integer(i4b),intent(out)        :: niter                         ! number of iterations taken
  integer(i4b),intent(out)        :: err                           ! error code
  character(*),intent(out)        :: message                       ! error message
  ! *********************************************************************************************************************************************************
  ! *********************************************************************************************************************************************************
  ! ---------------------------------------------------------------------------------------
  ! * general local variables
  ! ---------------------------------------------------------------------------------------
  character(LEN=256)              :: cmessage                      ! error message of downwind routine
  integer(i4b)                    :: iter                          ! iteration index
  integer(i4b)                    :: iVar                          ! index of variable
  integer(i4b)                    :: iLayer                        ! index of layer in the snow+soil domain
  integer(i4b)                    :: iState                        ! index of model state
  integer(i4b)                    :: nLeadDim                      ! length of the leading dimension of the Jacobian matrix (nBands or nState)
  integer(i4b)                    :: local_ixGroundwater           ! local index for groundwater representation
  real(rkind)                     :: bulkDensity                   ! bulk density of a given layer (kg m-3)
  real(rkind)                     :: volEnthalpy                   ! volumetric enthalpy of a given layer (J m-3)
  real(rkind),parameter           :: tempAccelerate=0.00_rkind        ! factor to force initial canopy temperatures to be close to air temperature
  real(rkind),parameter           :: xMinCanopyWater=0.0001_rkind     ! minimum value to initialize canopy water (kg m-2)
  real(rkind),parameter           :: tinyStep=0.000001_rkind          ! stupidly small time step (s)
  ! ------------------------------------------------------------------------------------------------------
  ! * model solver
  ! ------------------------------------------------------------------------------------------------------
  logical(lgt),parameter          :: forceFullMatrix=.true.       ! flag to force the use of the full Jacobian matrix
  integer(i4b)                    :: maxiter                       ! maximum number of iterations
  integer(i4b)                    :: ixMatrix                      ! form of matrix (band diagonal or full matrix)
  integer(i4b)                    :: localMaxIter                  ! maximum number of iterations (depends on solution type)
  integer(i4b), parameter         :: scalarMaxIter=100             ! maximum number of iterations for the scalar solution
  type(var_dlength)               :: flux_init                     ! model fluxes at the start of the time step
  real(rkind),allocatable         :: dBaseflow_dMatric(:,:)        ! derivative in baseflow w.r.t. matric head (s-1)  ! NOTE: allocatable, since not always needed
  real(rkind)                     :: stateVecNew(nState)           ! new state vector (mixed units)
  real(rkind)                     :: fluxVec0(nState)              ! flux vector (mixed units)
  real(rkind)                     :: fScale(nState)                ! characteristic scale of the function evaluations (mixed units)
  real(rkind)                     :: xScale(nState)                ! characteristic scale of the state vector (mixed units)
  real(rkind)                     :: dMat(nState)                  ! diagonal matrix (excludes flux derivatives)
  real(rkind)                     :: sMul(nState)    ! NOTE: qp    ! multiplier for state vector for the residual calculations
  real(rkind)                     :: rVec(nState)    ! NOTE: qp    ! residual vector
  real(rkind)                     :: rAdd(nState)                  ! additional terms in the residual vector
  real(rkind)                     :: fOld,fNew                     ! function values (-); NOTE: dimensionless because scaled
  real(rkind)                     :: xMin,xMax                     ! state minimum and maximum (mixed units)
  logical(lgt)                    :: converged                     ! convergence flag
  logical(lgt)                    :: feasible                      ! feasibility flag
  real(rkind)                     :: resSinkNew(nState)            ! additional terms in the residual vector
  real(rkind)                     :: fluxVecNew(nState)            ! new flux vector
  real(rkind)                     :: resVecNew(nState)  ! NOTE: qp ! new residual vector

  ! ---------------------------------------------------------------------------------------
  ! * sundials solver variables
  ! ---------------------------------------------------------------------------------------
  type(c_ptr)                     :: sunctx               ! SUNDIALS simulation context
  type(c_ptr)                     :: package_mem          ! SUNDIALS memory pointer 
  type(N_Vector),pointer          :: sunvec_y             ! SUNDIALS state vector
  type(N_Vector),pointer          :: sunvec_fscale        ! vector containing diagonal elements of scaling matrix
  type(N_Vector),pointer          :: sunvec_xscale        ! vector containing diagonal elements of scaling matrix
  type(SUNMatrix),pointer         :: sunmat_A             ! SUNDIALS Jacobian matrix
  type(SUNLinearSolver),pointer   :: sunlinsol_LS         ! sundials linear solver
  integer(i4b)                    :: retval, retvalr      ! return value
  integer(c_long)                 :: clong_nState         ! number of state variables but as a c_long  

  type(kinsol_data),target        :: kinsol_user_data     ! user data for the KINSOL solver
  real(c_double)                  :: fnormtol, scsteptol
  integer(c_long)                 :: mset
  real(c_double),dimension(nState):: scale
  ! ---------------------------------------------------------------------------------------
  ! point to variables in the data structures
  ! ---------------------------------------------------------------------------------------
  globalVars: associate(&
  ! model decisions
  ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,& ! intent(in):    [i4b]    groundwater parameterization
  ixSpatialGroundwater    => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,& ! intent(in):    [i4b]    spatial representation of groundwater (local-column or single-basin)
  ! check the need to merge snow layers
  mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in):    [dp(:)]  temperature of each snow/soil layer (K)
  mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in):    [dp(:)]  volumetric fraction of liquid water (-)
  mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(in):    [dp(:)]  volumetric fraction of ice (-)
  mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):    [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
  snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)         ,& ! intent(in):    [dp]     scaling parameter for the snow freezing curve (K-1)
  ! accelerate solution for temperature
  airtemp                 => forc_data%var(iLookFORCE%airtemp)                      ,& ! intent(in):    [dp]     temperature of the upper boundary of the snow and soil domains (K)
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy energy state variable
  ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy hydrology state variable (mass)
  ! vector of energy and hydrology indices for the snow and soil domains
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
  ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the soil domain
  nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):    [i4b]    number of energy state variables in the snow+soil domain
  nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):    [i4b]    number of hydrology state variables in the snow+soil domain
  nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in):    [i4b]    number of hydrology state variables in the soil domain
  ! mapping from full domain to the sub-domain
  ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):    [i4b]    mapping of full state vector to the state subset
  ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):    [i4b]    index of control volume for different domains (veg, snow, soil)
  ! type of state and domain for a given variable
  ixStateType_subset      => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in):    [i4b(:)] [state subset] type of desired model state variables
  ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat      ,& ! intent(in):    [i4b(:)] [state subset] domain for desired model state variables
  ! layer geometry
  nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):    [i4b]    number of snow layers
  nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):    [i4b]    number of soil layers
  nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in):    [i4b]    total number of layers
  )
  ! ---------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="systemSolv/"

  ! *****
  ! (0) PRELIMINARIES...
  ! ********************

  ! -----
  ! * initialize...
  ! ---------------

  ! check
  if(dt < tinyStep)then
    message=trim(message)//'dt is tiny'
    err=20; return
  endif

  ! initialize the flags
  tooMuchMelt        = .false.   ! too much melt
  reduceCoupledStep  = .false.   ! need to reduce the length of the coupled step

  ! define maximum number of iterations
  maxiter = nint(mpar_data%var(iLookPARAM%maxiter)%dat(1))

  ! modify the groundwater representation for this single-column implementation
  select case(ixSpatialGroundwater)
    case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
    case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
    case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
  end select ! (modify the groundwater representation for this single-column implementation)

  ! allocate space for the model fluxes at the start of the time step
  call allocLocal(flux_meta(:),flux_init,nSnow,nSoil,err,cmessage)
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  ! allocate space for the baseflow derivatives
  ! NOTE: needs allocation because only used when baseflow sinks are active
  if(ixGroundwater==qbaseTopmodel)then
    allocate(dBaseflow_dMatric(nSoil,nSoil),stat=err)  ! baseflow depends on total storage in the soil column, hence on matric head in every soil layer
    if(err/=0) print*, 'unable to allocate space for the baseflow derivatives non kinsol'
    allocate(kinsol_user_data%dBaseflow_dMatric(nSoil,nSoil),stat=err)
    if(err/=0) print*, 'unable to allocate space for the baseflow derivatives kinsol'
  else
    allocate(dBaseflow_dMatric(0,0),stat=err)          ! allocate zero-length dimnensions to avoid passing around an unallocated matrix
    if(err/=0) print*, 'unable to allocate space for the baseflow derivatives non kinsol'
    allocate(kinsol_user_data%dBaseflow_dMatric(0,0),stat=err)
    if(err/=0) print*, 'unable to allocate space for the baseflow derivatives non kinsol'
  end if
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the baseflow derivatives'; return; end if

  ! identify the matrix solution method
  ! (the type of matrix used to solve the linear system A.X=B)
  if(local_ixGroundwater==qbaseTopmodel .or. scalarSolution .or. forceFullMatrix)then
    nLeadDim=nState         ! length of the leading dimension
    ixMatrix=ixFullMatrix   ! named variable to denote the full Jacobian matrix
  else
    nLeadDim=nBands         ! length of the leading dimension
    ixMatrix=ixBandMatrix   ! named variable to denote the band-diagonal matrix
  endif
  ! initialize the model fluxes (some model fluxes are not computed in the iterations)
  do iVar=1,size(flux_temp%var)
    flux_init%var(iVar)%dat(:) = flux_temp%var(iVar)%dat(:)
  end do

    ! * get scaling vectors...
  ! ------------------------

  ! initialize state vectors
  call getScaling(&
                  ! input
                  diag_data,                        & ! intent(in):    model diagnostic variables for a local HRU
                  indx_data,                        & ! intent(in):    indices defining model states and layers
                  ! output
                  fScale,                           & ! intent(out):   function scaling vector (mixed units)
                  xScale,                           & ! intent(out):   variable scaling vector (mixed units)
                  sMul,                             & ! intent(out):   multiplier for state vector (used in the residual calculations)
                  dMat,                             & ! intent(out):   diagonal of the Jacobian matrix (excludes fluxes)
                  err,cmessage)                       ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! ***************************
  ! Set up the user data for kinsol
  ! ***************************
  kinsol_user_data%dt = dt
  kinsol_user_data%nSnow = nSnow 
  kinsol_user_data%nSoil = nSoil 
  kinsol_user_data%nLayers = nLayers 
  kinsol_user_data%computeVegFlux =  computeVegFlux
  kinsol_user_data%computeBaseflow = (ixGroundwater==qbaseTopmodel) 
  kinsol_user_data%ixMatrix = ixMatrix
  kinsol_user_data%nState = nState
  kinsol_user_data%scalarSolution = scalarSolution
  kinsol_user_data%firstSubStep = firstSubStep
  kinsol_user_data%firstFluxCall = firstFluxCall
  kinsol_user_data%firstSplitOper = firstSplitOper
  kinsol_user_data%ixSaturation = ixSaturation
  kinsol_user_data%feasible = feasible

  kinsol_user_data%type_data = type_data
  kinsol_user_data%attr_data = attr_data
  kinsol_user_data%mpar_data = mpar_data
  kinsol_user_data%forc_data = forc_data
  kinsol_user_data%bvar_data = bvar_data
  kinsol_user_data%flux_data = flux_temp

  kinsol_user_data%indx_data = indx_data 
  kinsol_user_data%prog_data = prog_data
  kinsol_user_data%diag_data = diag_data
  kinsol_user_data%deriv_data = deriv_data

  kinsol_user_data%dBaseflow_dMatric = dBaseflow_dMatric
  allocate(kinsol_user_data%dMat(nState), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the kinsol user data'; return; end if
  kinsol_user_data%dMat = dMat
  
  allocate(kinsol_user_data%sMul(nState), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the kinsol user data'; return; end if
  kinsol_user_data%sMul = sMul

  allocate(kinsol_user_data%fScale(nState), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the kinsol user data'; return; end if
  kinsol_user_data%fScale = fScale

  allocate(kinsol_user_data%fluxVec(nState), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the kinsol user data'; return; end if
  kinsol_user_data%fluxVec = fluxVec0

  allocate(kinsol_user_data%resSink(nState), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the kinsol user data'; return; end if
  kinsol_user_data%resSink = resSinkNew

  allocate(kinsol_user_data%model_decisions(maxvarDecisions), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the kinsol user data'; return; end if
  kinsol_user_data%model_decisions = model_decisions
  
  ! -----
  ! * compute the initial function evaluation...
  ! --------------------------------------------

  ! initialize the trial state vectors
  stateVecTrial = stateVecInit

  ! ***************************
  ! SUNDIALS Context Creation
  ! ***************************
  retval = FSUNContext_Create(c_null_ptr, sunctx)
  if(retval /= 0)then;err=20; message=trim(message)//'unable to create the SUNDIALS context';print*,message; return; endif



  ! Set the inital guess
  clong_nState = nState
  sunvec_y => FN_VMake_Serial(clong_nState, stateVecTrial, sunctx)
  if (.not. associated(sunvec_y)) then; err=20; message='systemSolv: sunvec = NULL'; print*,message;return; endif

  scale = 1.d0

  sunvec_fscale => FN_VMake_Serial(clong_nState, scale, sunctx)
  if (.not. associated(sunvec_fscale)) then; err=20; message='systemSolv: sunvec = NULL'; print*,message; return; endif

  package_mem = FKinCreate(sunctx)
  if (.not. c_associated(package_mem)) then; err=20; message='systemSolv: package_mem = NULL';print*,message; return; endif

  ! sunvec_xscale => FN_VMake_Serial(clong_nState, xScale, sunctx)
  ! if (.not. associated(sunvec_xscale)) then; err=20; message='systemSolv: sunvec = NULL'; return; endif

  ! Set the user data for sundials
  retval = FKINSetUserData(package_mem, c_loc(kinsol_user_data))
  if(retval /= 0)then; err=20; message=trim(message)//'unable to set the user data'; print*,message; return; endif

  retval = FKinInit(package_mem, c_funloc(eval8summa_kinsol), sunvec_y)
  if(retval /= 0)then; err=20; message=trim(message)//'unable to initialize the kinsol package'; print*,message; return; endif


  ! -------------------------
  ! Set optional inputs
  fnormtol = ftol
  retval = FKINSetFuncNormTol(package_mem, fnormtol)
  if(retval /= 0)then; err=20; message=trim(message)//'unable to set the function norm tolerance'; return; endif

  scsteptol = stol
  retval = FKINSetScaledStepTol(package_mem, scsteptol)
  if(retval /= 0)then; err=20; message=trim(message)//'unable to set the scaled step tolerance'; return; endif


  ! Create the Matrix Object
  ! Ax = b => A is the matrix, x is the state vector, b is the residual vector
  ! sunmat_A is the matrix
  sunmat_A => FSUNDenseMatrix(clong_nState, clong_nState, sunctx)
  if (.not. associated(sunmat_A)) then; err=20; message='systemSolv: sunmat = NULL'; return; endif

  ! Create the Linear Solver Object
  sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
  if (.not. associated(sunlinsol_LS)) then; err=20; message='systemSolv: sunlinsol = NULL'; return; endif

  ! Attach the matrix and linear solver to KINSOL
  retval = FKINSetLinearSolver(package_mem, sunlinsol_LS, sunmat_A)
  if(retval /= 0)then; err=20; message=trim(message)//'unable to set the linear solver'; return; endif

  ! Set up the Jacobian function
  retval = FKinSetJacFn(package_mem, c_funloc(computJacob_kinsol))
  if(retval /= 0)then; err=20; message=trim(message)//'unable to set the Jacobian function'; return; endif
  
  ! Indicate exact Newton
  mset = 1
  retval = FKINSetMaxSetupCalls(package_mem, mset)
  if(retval /= 0)then; err=20; message=trim(message)//'unable to set the maximum number of setup calls'; return; endif
  
  ! Call KINSol to solve problem
  retval = FKINSol(package_mem, sunvec_y, KIN_LINESEARCH, sunvec_fscale, sunvec_fscale)
  call free_kinsol_objects(package_mem, sunlinsol_LS, sunmat_A, sunvec_y, sunvec_fscale, sunctx)
  if(retval /= 0)then; err=20; message=trim(message)//'unable to solve the system of equations'; return; endif


  ! -----
  ! * update states...
  ! ------------------

  ! set untapped melt energy to zero
  untappedMelt(:) = 0._rkind

  ! update temperatures (ensure new temperature is consistent with the fluxes)
  if(nSnowSoilNrg>0)then
    do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
      iState = ixSnowSoilNrg(iLayer)
      stateVecTrial(iState) = stateVecInit(iState) + (fluxVecNew(iState)*dt + resSinkNew(iState))/real(sMul(iState), rkind)
    end do  ! looping through non-missing energy state variables in the snow+soil domain
  endif

  ! update volumetric water content in the snow (ensure change in state is consistent with the fluxes)
  ! NOTE: for soil water balance is constrained within the iteration loop
  if(nSnowSoilHyd>0)then
    do concurrent (iLayer=1:nSnow,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing water state variables in the snow domain)
      iState = ixSnowSoilHyd(iLayer)
      stateVecTrial(iState) = stateVecInit(iState) + (fluxVecNew(iState)*dt + resSinkNew(iState))
    end do  ! looping through non-missing water state variables in the soil domain
  endif

  ! end associate statements
  end associate globalVars

end subroutine systemSolv

subroutine free_kinsol_objects(package_mem, sunlinsol_LS, sunmat_A, sunvec_y, sunvec_fscale, sunctx)
  USE fsundials_context_mod                            ! Fortran interface to SUNContext
  USE fkinsol_mod                                      ! Fortran interface to KINSOL
  USE fsundials_nvector_mod                            ! Fortran interface to SUNDIALS N_Vector
  USE fsundials_matrix_mod                             ! Fortran interface to SUNDIALS Matrix
  USE fsundials_linearsolver_mod                       ! Fortran interface to generic SUNLinearSolver
  implicit none
  ! dummy variables
  type(c_ptr),           intent(inout) :: package_mem
  type(SUNLinearSolver), intent(inout) :: sunlinsol_LS
  type(SUNMatrix),       intent(inout) :: sunmat_A
  type(N_Vector),        intent(inout) :: sunvec_y
  type(N_Vector),        intent(inout) :: sunvec_fscale
  type(c_ptr),           intent(inout) :: sunctx
  ! local variables
  integer                              :: retval
  integer                              :: err
  character(len=256)                   :: message

  call FKINFree(package_mem)
  retval = FSUNLinSolFree(sunlinsol_LS)
  if(retval /= 0)then; err=20; message=trim(message)//'unable to free the linear solver'; print*,message; return; endif
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  ! call FN_VDestroy(sunvec_xscale)
  call FN_VDestroy(sunvec_fscale)
  retval = FSUNContext_Free(sunctx)
  if(retval /= 0)then; err=20; message=trim(message)//'unable to free the SUNDIALS context'; print*,message; return; endif
endsubroutine free_kinsol_objects

end module systemSolv_module
