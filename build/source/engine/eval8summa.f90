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

module eval8summa_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access the global print flag
USE globalData,only:globalPrintFlag

! define access to state variables to print
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! indices that define elements of the data structures
USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM                   ! named variables for structure elements
USE var_lookup,only:iLookPROG                    ! named variables for structure elements
USE var_lookup,only:iLookINDEX                   ! named variables for structure elements
USE var_lookup,only:iLookDIAG                    ! named variables for structure elements
USE var_lookup,only:iLookFLUX                    ! named variables for structure elements
USE var_lookup,only:iLookDERIV                   ! named variables for structure elements

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:  &
  localColumn,                & ! separate groundwater representation in each local soil column
  singleBasin                   ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  &
  qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
  bigBucket,                  & ! a big bucket (lumped aquifer model)
  noExplicit                    ! no explicit groundwater parameterization

! look-up values for the form of Richards' equation
USE mDecisions_module,only:  &
  moisture,                   & ! moisture-based form of Richards' equation
  mixdform                      ! mixed form of Richards' equation

implicit none
private
public::eval8summa
public::eval8summa_kinsol
public::write_residual_vector




contains

! **********************************************************************************************************
! public subroutine eval8summa: compute the residual vector and the Jacobian matrix
! **********************************************************************************************************
subroutine eval8summa(&
                      ! input: model control
                      dt,                      & ! intent(in):    length of the time step (seconds)
                      nSnow,                   & ! intent(in):    number of snow layers
                      nSoil,                   & ! intent(in):    number of soil layers
                      nLayers,                 & ! intent(in):    total number of layers
                      nState,                  & ! intent(in):    total number of state variables
                      firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                      firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                      firstSplitOper,          & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                      computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                      ! input: state vectors
                      stateVecTrial,           & ! intent(in):    model state vector
                      fScale,                  & ! intent(in):    function scaling vector
                      sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                      ! input: data structures
                      model_decisions,         & ! intent(in):    model decisions
                      type_data,               & ! intent(in):    type of vegetation and soil
                      attr_data,               & ! intent(in):    spatial attributes
                      mpar_data,               & ! intent(in):    model parameters
                      forc_data,               & ! intent(in):    model forcing data
                      bvar_data,               & ! intent(in):    average model variables for the entire basin
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      ! input-output: data structures
                      indx_data,               & ! intent(inout): index data
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,               & ! intent(inout): model fluxes for a local HRU
                      deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      ! input-output: baseflow
                      ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                      ! output: flux and residual vectors
                      feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                      fluxVec,                 & ! intent(out):   flux vector
                      resSink,                 & ! intent(out):   additional (sink) terms on the RHS of the state equation
                      resVec,                  & ! intent(out):   residual vector
                      fEval,                   & ! intent(out):   function evaluation
                      err,message)               ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! provide access to subroutines
  USE getVectorz_module, only:varExtract           ! extract variables from the state vector
  USE updateVars_module, only:updateVars           ! update prognostic variables
  USE computFlux_module, only:soilCmpres           ! compute soil compression
  USE computFlux_module, only:computFlux           ! compute fluxes given a state vector
  USE computResid_module,only:computResid          ! compute residuals given a state vector
  implicit none
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)             :: dt                     ! length of the time step (seconds)
  integer(i4b),intent(in)         :: nSnow                  ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                  ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                ! total number of layers
  integer(i4b),intent(in)         :: nState                 ! total number of state variables
  logical(lgt),intent(in)         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(inout)      :: firstFluxCall          ! flag to indicate if we are processing the first flux call
  logical(lgt),intent(in)         :: firstSplitOper         ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
  ! input: state vectors
  real(rkind),intent(in)             :: stateVecTrial(:)       ! model state vector
  real(rkind),intent(in)             :: fScale(:)              ! function scaling vector
  real(rkind),intent(in)             :: sMul(:)   ! NOTE: qp   ! state vector multiplier (used in the residual calculations)
  ! input: data structures
  type(model_options),intent(in)  :: model_decisions(:)     ! model decisions
  type(var_i),        intent(in)  :: type_data              ! type of vegetation and soil
  type(var_d),        intent(in)  :: attr_data              ! spatial attributes
  type(var_dlength),  intent(in)  :: mpar_data              ! model parameters
  type(var_d),        intent(in)  :: forc_data              ! model forcing data
  type(var_dlength),  intent(in)  :: bvar_data              ! model variables for the local basin
  type(var_dlength),  intent(in)  :: prog_data              ! prognostic variables for a local HRU
  ! output: data structures
  type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
  type(var_dlength),intent(inout) :: diag_data              ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_data              ! model fluxes for a local HRU
  type(var_dlength),intent(inout) :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
  ! input-output: baseflow
  integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
  real(rkind),intent(out)            :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
  ! output: flux and residual vectors
  logical(lgt),intent(out)        :: feasible               ! flag to denote the feasibility of the solution
  real(rkind),intent(out)            :: fluxVec(:)             ! flux vector
  real(rkind),intent(out)            :: resSink(:)             ! sink terms on the RHS of the flux equation
  real(rkind),intent(out)            :: resVec(:) ! NOTE: qp   ! residual vector
  real(rkind),intent(out)            :: fEval                  ! function evaluation
  ! output: error control
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! state variables
  real(rkind)                        :: scalarCanairTempTrial     ! trial value for temperature of the canopy air space (K)
  real(rkind)                        :: scalarCanopyTempTrial     ! trial value for temperature of the vegetation canopy (K)
  real(rkind)                        :: scalarCanopyWatTrial      ! trial value for liquid water storage in the canopy (kg m-2)
  real(rkind),dimension(nLayers)     :: mLayerTempTrial           ! trial value for temperature of layers in the snow and soil domains (K)
  real(rkind),dimension(nLayers)     :: mLayerVolFracWatTrial     ! trial value for volumetric fraction of total water (-)
  real(rkind),dimension(nSoil)       :: mLayerMatricHeadTrial     ! trial value for total water matric potential (m)
  real(rkind),dimension(nSoil)       :: mLayerMatricHeadLiqTrial  ! trial value for liquid water matric potential (m)
  real(rkind)                        :: scalarAquiferStorageTrial ! trial value of storage of water in the aquifer (m)
  ! diagnostic variables
  real(rkind)                        :: scalarCanopyLiqTrial      ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(rkind)                        :: scalarCanopyIceTrial      ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),dimension(nLayers)     :: mLayerVolFracLiqTrial     ! trial value for volumetric fraction of liquid water (-)
  real(rkind),dimension(nLayers)     :: mLayerVolFracIceTrial     ! trial value for volumetric fraction of ice (-)
  ! other local variables
  integer(i4b)                    :: iLayer                    ! index of model layer in the snow+soil domain
  integer(i4b)                    :: jState(1)                 ! index of model state for the scalar solution within the soil domain
  integer(i4b)                    :: ixBeg,ixEnd               ! index of indices for the soil compression routine
  integer(i4b),parameter          :: ixVegVolume=1             ! index of the desired vegetation control volumne (currently only one veg layer)
  real(rkind)                        :: xMin,xMax                 ! minimum and maximum values for water content
  real(rkind)                        :: scalarCanopyHydTrial      ! trial value for mass of water on the vegetation canopy (kg m-2)
  real(rkind),parameter              :: canopyTempMax=500._rkind     ! expected maximum value for the canopy temperature (K)
  real(rkind),dimension(nLayers)     :: mLayerVolFracHydTrial     ! trial value for volumetric fraction of water (-), general vector merged from Wat and Liq
  real(rkind),dimension(nState)      :: rVecScaled                ! scaled residual vector
  character(LEN=256)              :: cmessage                  ! error message of downwind routine
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! association to variables in the data structures
  ! --------------------------------------------------------------------------------------------------------------------------------
  associate(&
  ! model decisions
  ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision   ,&  ! intent(in):  [i4b]   index of the form of Richards' equation
  ! snow parameters
  snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)         ,&  ! intent(in):  [dp]    scaling parameter for the snow freezing curve (K-1)
  ! soil parameters
  theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,&  ! intent(in):  [dp(:)] soil porosity (-)
  specificStorage         => mpar_data%var(iLookPARAM%specificStorage)%dat(1)       ,&  ! intent(in):  [dp]    specific storage coefficient (m-1)
  ! canopy and layer depth
  canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,&  ! intent(in):  [dp   ] canopy depth (m)
  mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,&  ! intent(in):  [dp(:)] depth of each layer in the snow-soil sub-domain (m)
  ! model state variables
  scalarSfcMeltPond       => prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)      ,&  ! intent(in):  [dp]    ponded water caused by melt of the "snow without a layer" (kg m-2)
  mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,&  ! intent(in):  [dp(:)] volumetric fraction of liquid water (-)
  mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,&  ! intent(in):  [dp(:)] volumetric fraction of ice (-)
  mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,&  ! intent(in):  [dp(:)] liquid water matric potential (m)
  ! model diagnostic variables
  scalarFracLiqVeg        => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)       ,&  ! intent(in):  [dp]    fraction of liquid water on vegetation (-)
  mLayerFracLiqSnow       => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat         ,&  ! intent(in):  [dp(:)] fraction of liquid water in each snow layer (-)
  ! soil compression
  scalarSoilCompress      => diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1)     ,&  ! intent(in): [dp]    total change in storage associated with compression of the soil matrix (kg m-2)
  mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,&  ! intent(in): [dp(:)] change in storage associated with compression of the soil matrix (-)
  ! derivatives
  dVolTot_dPsi0           => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat           ,&  ! intent(in): [dp(:)] derivative in total water content w.r.t. total water matric potential
  dCompress_dPsi          => deriv_data%var(iLookDERIV%dCompress_dPsi)%dat          ,&  ! intent(in): [dp(:)] derivative in compressibility w.r.t. matric head (m-1)
  ! mapping
  ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,&  ! intent(in): [i4b(:)] mapping of full state vector to the state subset
  ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,&  ! intent(in): [i4b(:)] index of control volume for different domains (veg, snow, soil)
  ! indices
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,&  ! intent(in): [i4b]    index of canopy air space energy state variable (nrg)
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b]    index of canopy energy state variable (nrg)
  ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,&  ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
  ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
  ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,&  ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
  ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,&  ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
  ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,&  ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
  ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,&  ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                 &  ! intent(in): [i4b(:)] layer type (iname_soil or iname_snow)
  ) ! association to variables in the data structures
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="eval8summa/"

  ! check the feasibility of the solution
  feasible=.true.

  ! check that the canopy air space temperature is reasonable
  if(ixCasNrg/=integerMissing)then
  if(stateVecTrial(ixCasNrg) > canopyTempMax) feasible=.false.
  endif

  ! check that the canopy air space temperature is reasonable
  if(ixVegNrg/=integerMissing)then
  if(stateVecTrial(ixVegNrg) > canopyTempMax) feasible=.false.
  endif

  ! check canopy liquid water is not negative
  if(ixVegHyd/=integerMissing)then
  if(stateVecTrial(ixVegHyd) < 0._rkind) feasible=.false.
  end if

  ! check snow temperature is below freezing
  if(count(ixSnowOnlyNrg/=integerMissing)>0)then
  if(any(stateVecTrial( pack(ixSnowOnlyNrg,ixSnowOnlyNrg/=integerMissing) ) > Tfreeze)) feasible=.false.
  endif

  ! loop through non-missing hydrology state variables in the snow+soil domain
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)

  ! check the minimum and maximum water constraints
  if(ixHydType(iLayer)==iname_watLayer .or. ixHydType(iLayer)==iname_liqLayer)then

    ! --> minimum
    if (layerType(iLayer) == iname_soil) then
    xMin = theta_sat(iLayer-nSnow)
    else
    xMin = 0._rkind
    endif

    ! --> maximum
    select case( layerType(iLayer) )
    case(iname_snow); xMax = merge(iden_ice,  1._rkind - mLayerVolFracIce(iLayer), ixHydType(iLayer)==iname_watLayer)
    case(iname_soil); xMax = merge(theta_sat(iLayer-nSnow), theta_sat(iLayer-nSnow) - mLayerVolFracIce(iLayer), ixHydType(iLayer)==iname_watLayer)
    end select

    ! --> check
    if(stateVecTrial( ixSnowSoilHyd(iLayer) ) < xMin .or. stateVecTrial( ixSnowSoilHyd(iLayer) ) > xMax) feasible=.false.
    !if(.not.feasible) write(*,'(a,1x,i4,1x,L1,1x,10(f20.10,1x))') 'iLayer, feasible, stateVecTrial( ixSnowSoilHyd(iLayer) ), xMin, xMax = ', iLayer, feasible, stateVecTrial( ixSnowSoilHyd(iLayer) ), xMin, xMax

  endif  ! if water states

  end do  ! loop through non-missing hydrology state variables in the snow+soil domain

  ! early return for non-feasible solutions
  if(.not.feasible)then
    fluxVec(:) = realMissing
    resVec(:)  = quadMissing
    fEval      = realMissing
    print*, "Infeasible State"
    return
  end if

  ! get the start and end indices for the soil compression calculations
  if(scalarSolution)then
    jState = pack(ixControlVolume, ixMapFull2Subset/=integerMissing)
    ixBeg  = jState(1)
    ixEnd  = jState(1)
    else
    ixBeg  = 1
    ixEnd  = nSoil
  endif

  ! extract variables from the model state vector
  call varExtract(&
                  ! input
                  stateVecTrial,            & ! intent(in):    model state vector (mixed units)
                  diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                  prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                  indx_data,                & ! intent(in):    indices defining model states and layers
                  ! output: variables for the vegetation canopy
                  scalarCanairTempTrial,    & ! intent(out):   trial value of canopy air temperature (K)
                  scalarCanopyTempTrial,    & ! intent(out):   trial value of canopy temperature (K)
                  scalarCanopyWatTrial,     & ! intent(out):   trial value of canopy total water (kg m-2)
                  scalarCanopyLiqTrial,     & ! intent(out):   trial value of canopy liquid water (kg m-2)
                  scalarCanopyIceTrial,     & ! intent(out):   trial value of canopy ice content (kg m-2)
                  ! output: variables for the snow-soil domain
                  mLayerTempTrial,          & ! intent(out):   trial vector of layer temperature (K)
                  mLayerVolFracWatTrial,    & ! intent(out):   trial vector of volumetric total water content (-)
                  mLayerVolFracLiqTrial,    & ! intent(out):   trial vector of volumetric liquid water content (-)
                  mLayerVolFracIceTrial,    & ! intent(out):   trial vector of volumetric ice water content (-)
                  mLayerMatricHeadTrial,    & ! intent(out):   trial vector of total water matric potential (m)
                  mLayerMatricHeadLiqTrial, & ! intent(out):   trial vector of liquid water matric potential (m)
                  ! output: variables for the aquifer
                  scalarAquiferStorageTrial,& ! intent(out):   trial value of storage of water in the aquifer (m)
                  ! output: error control
                  err,cmessage)               ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

  ! update diagnostic variables
  call updateVars(&
                  ! input
                  .false.,                                   & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                  mpar_data,                                 & ! intent(in):    model parameters for a local HRU
                  indx_data,                                 & ! intent(in):    indices defining model states and layers
                  prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                                 & ! intent(inout): model diagnostic variables for a local HRU
                  deriv_data,                                & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                  ! output: variables for the vegetation canopy
                  scalarCanopyTempTrial,                     & ! intent(inout): trial value of canopy temperature (K)
                  scalarCanopyWatTrial,                      & ! intent(inout): trial value of canopy total water (kg m-2)
                  scalarCanopyLiqTrial,                      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                  scalarCanopyIceTrial,                      & ! intent(inout): trial value of canopy ice content (kg m-2)
                  ! output: variables for the snow-soil domain
                  mLayerTempTrial,                           & ! intent(inout): trial vector of layer temperature (K)
                  mLayerVolFracWatTrial,                     & ! intent(inout): trial vector of volumetric total water content (-)
                  mLayerVolFracLiqTrial,                     & ! intent(inout): trial vector of volumetric liquid water content (-)
                  mLayerVolFracIceTrial,                     & ! intent(inout): trial vector of volumetric ice water content (-)
                  mLayerMatricHeadTrial,                     & ! intent(inout): trial vector of total water matric potential (m)
                  mLayerMatricHeadLiqTrial,                  & ! intent(inout): trial vector of liquid water matric potential (m)
                  ! output: error control
                  err,cmessage)                                ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

  ! print the states in the canopy domain
  !print*, 'dt = ', dt
  !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyTempTrial    = ', scalarCanopyTempTrial
  !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyWatTrial     = ', scalarCanopyWatTrial
  !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyLiqTrial     = ', scalarCanopyLiqTrial
  !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyIceTrial     = ', scalarCanopyIceTrial

  ! print the states in the snow+soil domain
  !write(*,'(a,1x,10(f20.10,1x))') 'mLayerTempTrial          = ', mLayerTempTrial(iJac1:min(nLayers,iJac2))
  !write(*,'(a,1x,10(f20.10,1x))') 'mLayerVolFracWatTrial    = ', mLayerVolFracWatTrial(iJac1:min(nLayers,iJac2))
  !write(*,'(a,1x,10(f20.10,1x))') 'mLayerVolFracLiqTrial    = ', mLayerVolFracLiqTrial(iJac1:min(nLayers,iJac2))
  !write(*,'(a,1x,10(f20.10,1x))') 'mLayerVolFracIceTrial    = ', mLayerVolFracIceTrial(iJac1:min(nLayers,iJac2))
  !write(*,'(a,1x,10(f20.10,1x))') 'mLayerMatricHeadTrial    = ', mLayerMatricHeadTrial(iJac1:min(nSoil,iJac2))
  !write(*,'(a,1x,10(f20.10,1x))') 'mLayerMatricHeadLiqTrial = ', mLayerMatricHeadLiqTrial(iJac1:min(nSoil,iJac2))

  ! print the water content
  if(globalPrintFlag)then
  if(iJac1<nSnow) write(*,'(a,10(f16.10,1x))') 'mLayerVolFracWatTrial = ', mLayerVolFracWatTrial(iJac1:min(iJac2,nSnow))
  if(iJac1<nSnow) write(*,'(a,10(f16.10,1x))') 'mLayerVolFracLiqTrial = ', mLayerVolFracLiqTrial(iJac1:min(iJac2,nSnow))
  if(iJac1<nSnow) write(*,'(a,10(f16.10,1x))') 'mLayerVolFracIceTrial = ', mLayerVolFracIceTrial(iJac1:min(iJac2,nSnow))
  endif

  ! save the number of flux calls per time step
  indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) = indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) + 1

  ! compute the fluxes for a given state vector
  call computFlux(&
                  ! input-output: model control
                  nSnow,                     & ! intent(in):    number of snow layers
                  nSoil,                     & ! intent(in):    number of soil layers
                  nLayers,                   & ! intent(in):    total number of layers
                  firstSubStep,              & ! intent(in):    flag to indicate if we are processing the first sub-step
                  firstFluxCall,             & ! intent(inout): flag to denote the first flux call
                  firstSplitOper,            & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                  computeVegFlux,            & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                  scalarSolution,            & ! intent(in):    flag to indicate the scalar solution
                  scalarSfcMeltPond/dt,      & ! intent(in):    drainage from the surface melt pond (kg m-2 s-1)
                  ! input: state variables
                  scalarCanairTempTrial,     & ! intent(in):    trial value for the temperature of the canopy air space (K)
                  scalarCanopyTempTrial,     & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                  mLayerTempTrial,           & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                  mLayerMatricHeadLiqTrial,  & ! intent(in):    trial value for the liquid water matric potential in each soil layer (m)
                  scalarAquiferStorageTrial, & ! intent(in):    trial value of storage of water in the aquifer (m)
                  ! input: diagnostic variables defining the liquid water and ice content
                  scalarCanopyLiqTrial,      & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                  scalarCanopyIceTrial,      & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                  mLayerVolFracLiqTrial,     & ! intent(in):    trial value for the volumetric liquid water content in each snow and soil layer (-)
                  mLayerVolFracIceTrial,     & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                  ! input: data structures
                  model_decisions,           & ! intent(in):    model decisions
                  type_data,                 & ! intent(in):    type of vegetation and soil
                  attr_data,                 & ! intent(in):    spatial attributes
                  mpar_data,                 & ! intent(in):    model parameters
                  forc_data,                 & ! intent(in):    model forcing data
                  bvar_data,                 & ! intent(in):    average model variables for the entire basin
                  prog_data,                 & ! intent(in):    model prognostic variables for a local HRU
                  indx_data,                 & ! intent(in):    index data
                  ! input-output: data structures
                  diag_data,                 & ! intent(inout): model diagnostic variables for a local HRU
                  flux_data,                 & ! intent(inout): model fluxes for a local HRU
                  deriv_data,                & ! intent(out):   derivatives in model fluxes w.r.t. relevant state variables
                  ! input-output: flux vector and baseflow derivatives
                  ixSaturation,              & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                  dBaseflow_dMatric,         & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                  fluxVec,                   & ! intent(out):   flux vector (mixed units)
                  ! output: error control
                  err,cmessage)                ! intent(out):   error code and error message
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

  ! compute soil compressibility (-) and its derivative w.r.t. matric head (m)
  ! NOTE: we already extracted trial matrix head and volumetric liquid water as part of the flux calculations
  call soilCmpres(&
                  ! input:
                  ixRichards,                             & ! intent(in): choice of option for Richards' equation
                  ixBeg,ixEnd,                            & ! intent(in): start and end indices defining desired layers
                  mLayerMatricHeadLiq(1:nSoil),           & ! intent(in): matric head at the start of the time step (m)
                  mLayerMatricHeadLiqTrial(1:nSoil),      & ! intent(in): trial value of matric head (m)
                  mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in): trial value for the volumetric liquid water content in each soil layer (-)
                  mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in): trial value for the volumetric ice content in each soil layer (-)
                  specificStorage,                        & ! intent(in): specific storage coefficient (m-1)
                  theta_sat,                              & ! intent(in): soil porosity (-)
                  ! output:
                  mLayerCompress,                         & ! intent(inout): compressibility of the soil matrix (-)
                  dCompress_dPsi,                         & ! intent(inout): derivative in compressibility w.r.t. matric head (m-1)
                  err,cmessage)                             ! intent(out): error code and error message
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

  ! compute the total change in storage associated with compression of the soil matrix (kg m-2)
  scalarSoilCompress = sum(mLayerCompress(1:nSoil)*mLayerDepth(nSnow+1:nLayers))*iden_water

  ! vegetation domain: get the correct water states (total water, or liquid water, depending on the state type)
  if(computeVegFlux)then
  scalarCanopyHydTrial = merge(scalarCanopyWatTrial, scalarCanopyLiqTrial, (ixStateType( ixHydCanopy(ixVegVolume) )==iname_watCanopy) )
  else
  scalarCanopyHydTrial = realMissing
  endif

  ! snow+soil domain: get the correct water states (total water, or liquid water, depending on the state type)
  mLayerVolFracHydTrial = merge(mLayerVolFracWatTrial, mLayerVolFracLiqTrial, (ixHydType==iname_watLayer .or. ixHydType==iname_matLayer) )

  ! compute the residual vector
  call computResid(&
                  ! input: model control
                  dt,                        & ! intent(in):    length of the time step (seconds)
                  nSnow,                     & ! intent(in):    number of snow layers
                  nSoil,                     & ! intent(in):    number of soil layers
                  nLayers,                   & ! intent(in):    total number of layers
                  ! input: flux vectors
                  sMul,                      & ! intent(in):    state vector multiplier (used in the residual calculations)
                  fluxVec,                   & ! intent(in):    flux vector
                  ! input: state variables (already disaggregated into scalars and vectors)
                  scalarCanairTempTrial,     & ! intent(in):    trial value for the temperature of the canopy air space (K)
                  scalarCanopyTempTrial,     & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                  scalarCanopyHydTrial,      & ! intent(in):    trial value of canopy hydrology state variable (kg m-2)
                  mLayerTempTrial,           & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                  mLayerVolFracHydTrial,     & ! intent(in):    trial vector of volumetric water content (-)
                  scalarAquiferStorageTrial, & ! intent(in):    trial value of storage of water in the aquifer (m)
                  ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                  scalarCanopyIceTrial,      & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                  mLayerVolFracIceTrial,     & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                  ! input: data structures
                  prog_data,                 & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                 & ! intent(in):    model diagnostic variables for a local HRU
                  flux_data,                 & ! intent(in):    model fluxes for a local HRU
                  indx_data,                 & ! intent(in):    index data
                  ! output
                  resSink,                   & ! intent(out):   additional (sink) terms on the RHS of the state equation
                  resVec,                    & ! intent(out):   residual vector
                  err,cmessage)                ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

  ! compute the function evaluation
  rVecScaled = fScale(:)*real(resVec(:), rkind)   ! scale the residual vector (NOTE: residual vector is in quadruple precision)
  fEval      = 0.5_rkind*dot_product(rVecScaled,rVecScaled)

  call write_residual_vector(resVec, nState, message, err)

  ! end association with the information in the data structures
  end associate

end subroutine eval8summa



integer(c_int) function eval8summa_kinsol(sunvec_y, sunvec_f, user_data) &
                            result(ierr) bind(C,name='eval8summa_kinsol')
  USE, intrinsic :: iso_c_binding
  USE fsundials_nvector_mod
  USE fnvector_serial_mod
  USE kinsol_user_data_type

  implicit none
  type(N_Vector)             :: sunvec_y  ! solution N_Vector
  type(N_Vector)             :: sunvec_f  ! rhs N_Vector
  type(c_ptr),   value       :: user_data ! user-defined data

  ! local variables
  real(c_double),pointer     :: stateVec(:)
  real(c_double),pointer     :: resVec(:)
  type(kinsol_data), pointer :: kinsol_user_data
  integer(i4b)               :: err
  character(len=256)         :: message

  ierr = 0

  call c_f_pointer(user_data, kinsol_user_data)  ! associate the user data with the data structure
  stateVec(1:kinsol_user_data%nState) => FN_VGetArrayPointer(sunvec_y)
  resVec(1:kinsol_user_data%nState)   => FN_VGetArrayPointer(sunvec_f)


  if (kinsol_user_data%firstStateiteration) then
    kinsol_user_data%firstStateiteration = .false.
  else
    call imposeConstraints(kinsol_user_data%indx_data,kinsol_user_data%prog_data,kinsol_user_data%mpar_data,stateVec(:), &
       kinsol_user_data%stateVecprev, kinsol_user_data%nState, kinsol_user_data%nSoil, kinsol_user_data%nSnow, message, err)
    if(err/=0)then; ierr=1; message="Impose Constraints Failed"; print*, message; return; end if  ! (check for errors)

    kinsol_user_data%stateVecprev = stateVec(:)  
  endif

  call eval8summa(kinsol_user_data%dt,                &
                  kinsol_user_data%nSnow,             &
                  kinsol_user_data%nSoil,             &
                  kinsol_user_data%nLayers,           &
                  kinsol_user_data%nState,            &
                  kinsol_user_data%firstSubStep,      &
                  kinsol_user_data%firstFluxCall,     &
                  kinsol_user_data%firstSplitOper,    &
                  kinsol_user_data%computeVegFlux,    &
                  kinsol_user_data%scalarSolution,    &
                  stateVec,                           &
                  kinsol_user_data%fScale,            &
                  kinsol_user_data%sMul,              &
                  kinsol_user_data%model_decisions,   &
                  kinsol_user_data%type_data,         &
                  kinsol_user_data%attr_data,         &
                  kinsol_user_data%mpar_data,         &
                  kinsol_user_data%forc_data,         &
                  kinsol_user_data%bvar_data,         &
                  kinsol_user_data%prog_data,         &
                  kinsol_user_data%indx_data,         &
                  kinsol_user_data%diag_data,         &
                  kinsol_user_data%flux_data,         &
                  kinsol_user_data%deriv_data,        &
                  kinsol_user_data%ixSaturation,      &
                  kinsol_user_data%dBaseflow_dMatric, &
                  kinsol_user_data%feasible,          &
                  kinsol_user_data%fluxVec,           &
                  kinsol_user_data%resSink,           &
                  resVec,                             &
                  kinsol_user_data%fEval,             &
                  ierr, message)

end function eval8summa_kinsol

subroutine imposeConstraints(indx_data, prog_data, mpar_data, stateVec, stateVecPrev,&
    nState, nSoil, nSnow, message, err)
      ! external functions
  USE snow_utils_module,only:fracliquid                           ! compute the fraction of liquid water at a given temperature (snow)
  USE soil_utils_module,only:crit_soilT                           ! compute the critical temperature below which ice exists
  
  implicit none
  type(var_ilength),intent(in)             :: indx_data                   ! indices defining model states and layers
  type(var_dlength),intent(in)             :: prog_data                   ! prognostic variables for a local HRU
  type(var_dlength),intent(in)             :: mpar_data                   ! model parameters
  real(rkind), intent(inout)               :: stateVec(:)
  real(rkind), intent(in)                  :: stateVecPrev(:)
  integer(i4b), intent(in)                 :: nState
  integer(i4b), intent(in)                 :: nSoil
  integer(i4b), intent(in)                 :: nSnow
  character(len=256), intent(out)          :: message
  integer(i4b), intent(out)                :: err

  real(rkind),dimension(nState)            :: xInc                     ! iteration increment
  ! -----------------------------------------------------------------------------------------------------
  ! temporary variables for model constraints
  real(rkind)                              :: cInc                         ! constrained temperature increment (K) -- simplified bi-section
  real(rkind)                              :: xIncFactor                   ! scaling factor for the iteration increment (-)
  integer(i4b)                             :: iMax(1)                      ! index of maximum temperature
  real(rkind)                              :: scalarTemp                   ! temperature of an individual snow layer (K)
  real(rkind)                              :: volFracLiq                   ! volumetric liquid water content of an individual snow layer (-)
  logical(lgt),dimension(nSoil)            :: crosFlag                     ! flag to denote temperature crossing from unfrozen to frozen (or vice-versa)
  logical(lgt)                             :: crosTempVeg                  ! flag to denoote where temperature crosses the freezing point
  real(rkind)                              :: xPsi00                       ! matric head after applying the iteration increment (m)
  real(rkind)                              :: TcSoil                       ! critical point when soil begins to freeze (K)
  real(rkind)                              :: critDiff                     ! temperature difference from critical (K)
  real(rkind),parameter                    :: epsT=1.e-7_rkind                ! small interval above/below critical (K)
  real(rkind),parameter                    :: zMaxTempIncrement=1._rkind      ! maximum temperature increment (K)
  ! indices of model state variables
  integer(i4b)                             :: iState                       ! index of state within a specific variable type
  integer(i4b)                             :: ixNrg,ixLiq                  ! index of energy and mass state variables in full state vector
  ! indices of model layers
  integer(i4b)                             :: iLayer                       ! index of model layer
  ! -----------------------------------------------------------------------------------------------------
  ! associate variables with indices of model state variables
  associate(&
  ixNrgOnly               => indx_data%var(iLookINDEX%ixNrgOnly)%dat                ,& ! intent(in): [i4b(:)] list of indices in the state subset for energy states
  ixHydOnly               => indx_data%var(iLookINDEX%ixHydOnly)%dat                ,& ! intent(in): [i4b(:)] list of indices in the state subset for hydrology states
  ixMatOnly               => indx_data%var(iLookINDEX%ixMatOnly)%dat                ,& ! intent(in): [i4b(:)] list of indices in the state subset for matric head states
  ixMassOnly              => indx_data%var(iLookINDEX%ixMassOnly)%dat               ,& ! intent(in): [i4b(:)] list of indices in the state subset for canopy storage states
  ixStateType_subset      => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in): [i4b(:)] named variables defining the states in the subset
  ! indices for specific state variables
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in): [i4b] index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in): [i4b] index of canopy energy state variable
  ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in): [i4b] index of canopy hydrology state variable (mass)
  ixTopNrg                => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)              ,& ! intent(in): [i4b] index of upper-most energy state in the snow-soil subdomain
  ixTopHyd                => indx_data%var(iLookINDEX%ixTopHyd)%dat(1)              ,& ! intent(in): [i4b] index of upper-most hydrology state in the snow-soil subdomain
  ! vector of energy indices for the snow and soil domains
  ! NOTE: states not in the subset are equal to integerMissing
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
  ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow domain
  ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the soil domain
  ! vector of hydrology indices for the snow and soil domains
  ! NOTE: states not in the subset are equal to integerMissing
  ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
  ixSnowOnlyHyd           => indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow domain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the soil domain
  ! number of state variables of a specific type
  nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
  nSnowOnlyNrg            => indx_data%var(iLookINDEX%nSnowOnlyNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the snow domain
  nSoilOnlyNrg            => indx_data%var(iLookINDEX%nSoilOnlyNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
  nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
  nSnowOnlyHyd            => indx_data%var(iLookINDEX%nSnowOnlyHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
  nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
  ! state variables at the start of the time step
  mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat           & ! intent(in): [dp(:)] matric head (m)
  ) ! associating variables with indices of model state variables
  ! -----------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='imposeConstraints/'

  xInc(:) = stateVec(:) - stateVecPrev(:)

  ! ** limit temperature increment to zMaxTempIncrement
  if(any(abs(xInc(ixNrgOnly)) > zMaxTempIncrement))then
    iMax       = maxloc( abs(xInc(ixNrgOnly)) )                     ! index of maximum temperature increment
    xIncFactor = abs( zMaxTempIncrement/xInc(ixNrgOnly(iMax(1))) )  ! scaling factor for the iteration increment (-)
    xInc       = xIncFactor*xInc
  end if

  ! ** impose solution constraints for vegetation
  ! (stop just above or just below the freezing point if crossing)
  ! --------------------------------------------------------------------------------------------------------------------
  ! canopy temperatures

  if(ixVegNrg/=integerMissing)then

    ! initialize
    critDiff    = Tfreeze - stateVecPrev(ixVegNrg)
    crosTempVeg = .false.

    ! initially frozen (T < Tfreeze)
    if(critDiff > 0._rkind)then
      if(xInc(ixVegNrg) > critDiff)then
        crosTempVeg = .true.
        cInc        = critDiff + epsT  ! constrained temperature increment (K)
      end if

    ! initially unfrozen (T > Tfreeze)
    else
      if(xInc(ixVegNrg) < critDiff)then
        crosTempVeg = .true.
        cInc        = critDiff - epsT  ! constrained temperature increment (K)
      end if

  end if  ! switch between frozen and unfrozen

  ! scale iterations
  if(crosTempVeg)then
    xIncFactor  = cInc/xInc(ixVegNrg)  ! scaling factor for the iteration increment (-)
    xInc        = xIncFactor*xInc      ! scale iteration increments
  endif

  endif  ! if the state variable for canopy temperature is included within the state subset

  ! --------------------------------------------------------------------------------------------------------------------
  ! canopy liquid water

  if(ixVegHyd/=integerMissing)then

    ! check if new value of storage will be negative
    if(stateVecPrev(ixVegHyd)+xInc(ixVegHyd) < 0._rkind)then
      ! scale iteration increment
      cInc       = -0.5_rkind*stateVecPrev(ixVegHyd)              ! constrained iteration increment (K) -- simplified bi-section
      xIncFactor = cInc/xInc(ixVegHyd)                             ! scaling factor for the iteration increment (-)
      xInc       = xIncFactor*xInc                                  ! new iteration increment
    end if

  endif  ! if the state variable for canopy water is included within the state subset

  ! --------------------------------------------------------------------------------------------------------------------
  ! ** impose solution constraints for snow
  if(nSnowOnlyNrg > 0)then

    ! loop through snow layers
    checksnow: do iLayer=1,nSnow  ! necessary to ensure that NO layers rise above Tfreeze

      ! check of the data is mising
      if(ixSnowOnlyNrg(iLayer)==integerMissing) cycle

      ! check temperatures, and, if necessary, scale iteration increment
      iState = ixSnowOnlyNrg(iLayer)
      if(stateVecPrev(iState) + xInc(iState) > Tfreeze)then
      ! scale iteration increment
        cInc       = 0.5_rkind*(Tfreeze - stateVecPrev(iState) )        ! constrained temperature increment (K) -- simplified bi-section
        xIncFactor = cInc/xInc(iState)                                ! scaling factor for the iteration increment (-)
        xInc       = xIncFactor*xInc
      end if   ! if snow temperature > freezing

    end do checkSnow

  endif  ! if there are state variables for energy in the snow domain

  ! --------------------------------------------------------------------------------------------------------------------
  ! - check if drain more than what is available
  ! NOTE: change in total water is only due to liquid flux
  if(nSnowOnlyHyd>0)then

    ! loop through snow layers
    do iLayer=1,nSnow

      ! * check if the layer is included
      if(ixSnowOnlyHyd(iLayer)==integerMissing) cycle

      ! * get the layer temperature (from stateVecPrev if ixSnowOnlyNrg(iLayer) is within the state vector
      if(ixSnowOnlyNrg(iLayer)/=integerMissing)then
        scalarTemp = stateVecPrev( ixSnowOnlyNrg(iLayer) )

      ! * get the layer temperature from the last update
      else
        scalarTemp = prog_data%var(iLookPROG%mLayerTemp)%dat(iLayer)
      endif

      ! * get the volumetric fraction of liquid water
      select case( ixStateType_subset( ixSnowOnlyHyd(iLayer) ) )
        case(iname_watLayer); volFracLiq = fracliquid(scalarTemp,mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)) * stateVecPrev(ixSnowOnlyHyd(iLayer))
        case(iname_liqLayer); volFracLiq = stateVecPrev(ixSnowOnlyHyd(iLayer))
        case default; err=20; message=trim(message)//'expect ixStateType_subset to be iname_watLayer or iname_liqLayer for snow hydrology'; return
      end select

      ! * check that the iteration increment does not exceed volumetric liquid water content
      if(-xInc(ixSnowOnlyHyd(iLayer)) > volFracLiq)then
        xInc(ixSnowOnlyHyd(iLayer)) = -0.5_rkind*volFracLiq
      endif

    end do  ! looping through snow layers

  endif   ! if there are state variables for liquid water in the snow domain

  ! --------------------------------------------------------------------------------------------------------------------
  ! ** impose solution constraints for soil temperature
  if(nSoilOnlyNrg>0)then
    do iLayer=1,nSoil

      ! - check if energy state is included
      if(ixSoilOnlyNrg(iLayer)==integerMissing) cycle

      ! - define index of the state variables within the state subset
      ixNrg = ixSoilOnlyNrg(iLayer)
      ixLiq = ixSoilOnlyHyd(iLayer)

      ! get the matric potential of total water
      if(ixLiq/=integerMissing)then
        xPsi00 = stateVecPrev(ixLiq) + xInc(ixLiq)
      else
        xPsi00 = mLayerMatricHead(iLayer)
      endif

      ! identify the critical point when soil begins to freeze (TcSoil)
      TcSoil = crit_soilT(xPsi00)

      ! get the difference from the current state and the crossing point (K)
      critDiff = TcSoil - stateVecPrev(ixNrg)

      ! * initially frozen (T < TcSoil)
      if(critDiff > 0._rkind)then

        ! (check crossing above zero)
        if(xInc(ixNrg) > critDiff)then
          ! print*, "xInc(ixNrg) > critDiff"
          crosFlag(iLayer) = .true.
          xInc(ixNrg) = critDiff + epsT  ! set iteration increment to slightly above critical temperature
        endif

      ! * initially unfrozen (T > TcSoil)
      else

        ! (check crossing below zero)
        if(xInc(ixNrg) < critDiff)then
          ! print*, "xInc(ixNrg) < critDiff"
          crosFlag(iLayer) = .true.
          xInc(ixNrg) = critDiff - epsT  ! set iteration increment to slightly below critical temperature
        endif

      endif  ! (switch between initially frozen and initially unfrozen)

    end do  ! (loop through soil layers)
  endif   ! (if there are both energy and liquid water state variables)

  ! ** impose solution constraints matric head
  if(size(ixMatOnly)>0)then
    do iState=1,size(ixMatOnly)

      ! - define index of the hydrology state variable within the state subset
      ixLiq = ixMatOnly(iState)

      ! - place constraint for matric head
      if(xInc(ixLiq) > 1._rkind .and. stateVecPrev(ixLiq) > 0._rkind)then
        xInc(ixLiq) = 1._rkind
      endif  ! if constraining matric head

    end do  ! (loop through soil layers)
  endif   ! (if there are both energy and liquid water state variables)

  ! end association with variables with indices of model state variables
  end associate

  ! Update the state vector with the modified iteration increment
  stateVec(:) = stateVecPrev(:) + xInc(:)


end subroutine imposeConstraints

subroutine write_residual_vector(resVec,nState,message,err)
  USE csv_file
  implicit none
  real(rkind),        intent(in)  :: resVec(:)
  integer(i4b),       intent(in)  :: nState
  character(len=256), intent(out) :: message
  integer(i4b),       intent(out) :: err

  logical(lgt)      :: fileExists       


  ! check if the file exists
  inquire(file='data_kinsol.csv', exist=fileExists)
  if (.not. fileExists) then
    open(unit=10, file ='data_kinsol.csv', status="new")
      call csv_write(10, resVec(1:nState), .true.)  
  else
    open(unit=10, file ='data_kinsol.csv', status="old", position="append")
      call csv_write(10, resVec(1:nState), .true.)  
  end if

  close(10)
endsubroutine write_residual_vector





end module eval8summa_module
