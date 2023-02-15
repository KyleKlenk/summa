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

module computJacob_module
USE, intrinsic :: iso_c_binding

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (dp)

! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements


! access missing values
USE globalData,only:integerMissing  ! missing integer


! access named variables to describe the form and structure of the matrices used in the numerical solver
USE globalData,only: ku             ! number of super-diagonal bands
USE globalData,only: kl             ! number of sub-diagonal bands
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

implicit none
! define constants
real(rkind),parameter     :: verySmall=tiny(1.0_rkind)     ! a very small number

private
public::computJacob
public::computJacob_kinsol
contains

 ! **********************************************************************************************************
 ! public subroutine computJacob: compute the Jacobian matrix
 ! **********************************************************************************************************
subroutine computJacob(&
                        ! input: model control
                        dt,                         & ! intent(in):    length of the time step (seconds)
                        nSnow,                      & ! intent(in):    number of snow layers
                        nSoil,                      & ! intent(in):    number of soil layers
                        nLayers,                    & ! intent(in):    total number of layers
                        computeVegFlux,             & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                        computeBaseflow,            & ! intent(in):    flag to indicate if we need to compute baseflow
                        ixMatrix,                   & ! intent(in):    form of the Jacobian matrix
                        ! input: data structures
                        indx_data,                  & ! intent(in):    index data
                        prog_data,                  & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,                  & ! intent(in):    model diagnostic variables for a local HRU
                        deriv_data,                 & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                        dBaseflow_dMatric,          & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                        ! input-output: Jacobian and its diagonal
                        dMat,                       & ! intent(inout): diagonal of the Jacobian matrix
                        aJac,                       & ! intent(out):   Jacobian matrix
                        ! output: error control
                        err,message)                  ! intent(out):   error code and error message
  ! -----------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  real(rkind),intent(in)            :: dt              ! length of the time step (seconds)
  integer(i4b),intent(in)           :: nSnow           ! number of snow layers
  integer(i4b),intent(in)           :: nSoil           ! number of soil layers
  integer(i4b),intent(in)           :: nLayers         ! total number of layers in the snow+soil domain
  logical(lgt),intent(in)           :: computeVegFlux  ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)           :: computeBaseflow ! flag to indicate if computing baseflow
  integer(i4b),intent(in)           :: ixMatrix        ! form of the Jacobian matrix
  ! input: data structures
  type(var_ilength),intent(in)      :: indx_data       ! indices defining model states and layers
  type(var_dlength),intent(in)      :: prog_data       ! prognostic variables for a local HRU
  type(var_dlength),intent(in)      :: diag_data       ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)      :: deriv_data      ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind),intent(in)            :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
  ! input-output: Jacobian and its diagonal
  real(rkind),intent(inout)         :: dMat(:)         ! diagonal of the Jacobian matrix
  real(rkind),intent(out)           :: aJac(:,:)       ! Jacobian matrix
  ! output variables
  integer(i4b),intent(out)          :: err             ! error code
  character(*),intent(out)          :: message         ! error message
  ! --------------------------------------------------------------
  ! * local variables
  ! --------------------------------------------------------------
  ! indices of model state variables
  integer(i4b)                      :: jState          ! index of state within the state subset
  integer(i4b)                      :: qState          ! index of cross-derivative state variable for baseflow
  integer(i4b)                      :: nrgState        ! energy state variable
  integer(i4b)                      :: watState        ! hydrology state variable
  integer(i4b)                      :: nState          ! number of state variables
  ! indices of model layers
  integer(i4b)                      :: iLayer          ! index of model layer
  integer(i4b)                      :: jLayer          ! index of model layer within the full state vector (hydrology)
  integer(i4b)                      :: pLayer          ! indices of soil layers (used for the baseflow derivatives)
  ! conversion factors
  real(rkind)                          :: convLiq2tot     ! factor to convert liquid water derivative to total water derivative
    ! --------------------------------------------------------------
  ! associate variables from data structures
  associate(&
  ! indices of model state variables
  ixTopNrg                     => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)                       ,& ! intent(in): [i4b] index of upper-most energy state in the snow+soil subdomain
  ! vectors of indices for specfic state types within specific sub-domains IN THE FULL STATE VECTOR
  ixNrgLayer                   => indx_data%var(iLookINDEX%ixNrgLayer)%dat                        ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
  ! vector of energy indices for the snow and soil domains
  ! NOTE: states not in the subset are equal to integerMissing
  ixSnowSoilNrg                => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
  ixSoilOnlyNrg                => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the soil domain
  ! vector of hydrology indices for the snow and soil domains
  ! NOTE: states not in the subset are equal to integerMissing
  ixSnowSoilHyd                => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
  ixSoilOnlyHyd                => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the soil domain
  ! number of state variables of a specific type
  nSnowSoilNrg                 => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
  nSoilOnlyNrg                 => indx_data%var(iLookINDEX%nSoilOnlyNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
  nSnowSoilHyd                 => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
  nSoilOnlyHyd                 => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
  ! derivatives in evaporative fluxes w.r.t. relevant state variables
  dGroundEvaporation_dTGround  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTGround )%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. ground temperature
  ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. temperature in layers above and below
  dNrgFlux_dTempAbove          => deriv_data%var(iLookDERIV%dNrgFlux_dTempAbove         )%dat     ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. temperature in the layer above
  dNrgFlux_dTempBelow          => deriv_data%var(iLookDERIV%dNrgFlux_dTempBelow         )%dat     ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. temperature in the layer below
  ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
  dVolTot_dPsi0                => deriv_data%var(iLookDERIV%dVolTot_dPsi0               )%dat     ,& ! intent(in): [dp(:)]  derivative in total water content w.r.t. total water matric potential
  dq_dHydStateAbove            => deriv_data%var(iLookDERIV%dq_dHydStateAbove           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer above
  dq_dHydStateBelow            => deriv_data%var(iLookDERIV%dq_dHydStateBelow           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer below
  dCompress_dPsi               => deriv_data%var(iLookDERIV%dCompress_dPsi              )%dat     ,& ! intent(in): [dp(:)]  derivative in compressibility w.r.t matric head
  ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
  dq_dNrgStateAbove            => deriv_data%var(iLookDERIV%dq_dNrgStateAbove           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer above
  dq_dNrgStateBelow            => deriv_data%var(iLookDERIV%dq_dNrgStateBelow           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer below
  mLayerdTheta_dTk             => deriv_data%var(iLookDERIV%mLayerdTheta_dTk            )%dat     ,& ! intent(in): [dp(:)]  derivative of volumetric liquid water content w.r.t. temperature
  ! diagnostic variables
  mLayerVolHtCapBulk           => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat                 ,& ! intent(in): [dp(:)]  bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
  ! layer depth
  mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth)%dat                         & ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
  ) ! making association with data in structures
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='computJacob/'

  ! *********************************************************************************************************************************************************
  ! *********************************************************************************************************************************************************
  ! * PART 0: PRELIMINARIES (INITIALIZE JACOBIAN AND COMPUTE TIME-VARIABLE DIAGONAL TERMS)
  ! *********************************************************************************************************************************************************
  ! *********************************************************************************************************************************************************

  ! get the number of state variables
  nState = size(dMat)

  ! initialize the Jacobian
  ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver
  aJac(:,:) = 0._rkind  ! analytical Jacobian matrix

  ! compute additional terms for the Jacobian for the snow-soil domain (excluding fluxes)
  ! NOTE: energy for snow+soil is computed *within* the iteration loop as it includes phase change
  do iLayer=1,nLayers
    if(ixSnowSoilNrg(iLayer)/=integerMissing) dMat(ixSnowSoilNrg(iLayer)) = mLayerVolHtCapBulk(iLayer) + LH_fus*iden_water*mLayerdTheta_dTk(iLayer)
  end do

  ! compute additional terms for the Jacobian for the soil domain (excluding fluxes)
  do iLayer=1,nSoil
    if(ixSoilOnlyHyd(iLayer)/=integerMissing) dMat(ixSoilOnlyHyd(iLayer)) = dVolTot_dPsi0(iLayer) + dCompress_dPsi(iLayer)
  end do

  ! define the form of the matrix
  select case(ixMatrix)
    case(ixBandMatrix)
      err = 10; message=trim(message)//'band matrix not implemented'; print*, message; return

    case(ixFullMatrix)
      ! check
      if(size(aJac,1)/=size(dMat) .or. size(aJac,2)/=size(dMat))then
        message=trim(message)//'unexpected shape of the Jacobian matrix: expect aJac(nState,nState)'
        err=20; return
      end if

      ! -----
      ! * energy fluxes for the snow+soil domain...
      ! -------------------------------------------
      if(nSnowSoilNrg>0)then
        do iLayer=1,nLayers  ! loop through all layers in the snow+soil domain

          ! check if the state is in the subset
          if(ixSnowSoilNrg(iLayer)==integerMissing) cycle
            
          ! - define index within the state subset and the full state vector
          jState = ixSnowSoilNrg(iLayer)        ! index within the state subset

          ! - diagonal elements
          aJac(jState,jState)   = (dt/mLayerDepth(iLayer))*(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + dMat(jState)

          ! - lower-diagonal elements
          if(iLayer > 1)then
            if(ixSnowSoilNrg(iLayer-1)/=integerMissing) aJac(ixSnowSoilNrg(iLayer-1),jState) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
          endif

          ! - upper diagonal elements
          if(iLayer < nLayers)then
            if(ixSnowSoilNrg(iLayer+1)/=integerMissing) aJac(ixSnowSoilNrg(iLayer+1),jState) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
          endif

        end do  ! (looping through energy states in the snow+soil domain)
      endif   ! (if the subset includes energy state variables in the snow+soil domain)

      ! -----
      ! * liquid water fluxes for the soil domain...
      ! --------------------------------------------
      if(nSoilOnlyHyd>0)then

        do iLayer=1,nSoil

          ! - check that the soil layer is desired
          if(ixSoilOnlyHyd(iLayer)==integerMissing) cycle

          ! - define state indices
          watState = ixSoilOnlyHyd(iLayer)         ! hydrology state index within the state subset

          ! - define indices of the soil layers
          jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector

          ! - compute the diagonal elements
          ! all terms *excluding* baseflow
          aJac(watState,watState) = (dt/mLayerDepth(jLayer))*(-dq_dHydStateBelow(iLayer-1) + dq_dHydStateAbove(iLayer)) + dMat(watState)

          ! - compute the lower-diagonal elements
          if(iLayer > 1)then
            if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer-1),watState) = (dt/mLayerDepth(jLayer-1))*( dq_dHydStateBelow(iLayer-1))
          endif

          ! - compute the upper-diagonal elements
          if(iLayer<nSoil)then
            if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer+1),watState) = (dt/mLayerDepth(jLayer+1))*(-dq_dHydStateAbove(iLayer))
          endif

        end do  ! (looping through hydrology states in the soil domain)
      endif   ! (if the subset includes hydrology state variables in the soil domain)

      ! -----
      ! * derivative in liquid water fluxes w.r.t. temperature for the soil domain...
      ! -----------------------------------------------------------------------------
      if(nSoilOnlyHyd>0 .and. nSoilOnlyNrg>0)then
        do iLayer=1,nSoilOnlyHyd

        ! - check that the soil layer is desired
        if(ixSoilOnlyHyd(iLayer)==integerMissing) cycle

        ! - define index of hydrology state variable within the state subset
        watState = ixSoilOnlyHyd(iLayer)

        ! - define indices of the soil layers
        jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector

        ! - define the energy state variable
        nrgState = ixNrgLayer(jLayer)       ! index within the full state vector

        ! only compute derivatives if the energy state for the current layer is within the state subset
        if(nrgstate/=integerMissing)then

          ! - compute the Jacobian for the layer itself
          aJac(watState,nrgState) = (dt/mLayerDepth(jLayer))*(-dq_dNrgStateBelow(iLayer-1) + dq_dNrgStateAbove(iLayer))   ! dVol/dT (K-1) -- flux depends on ice impedance

          ! - include derivatives w.r.t. ground evaporation
          if(nSnow==0 .and. iLayer==1)then  ! upper-most soil layer
            aJac(watState,ixTopNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTGround/iden_water) + aJac(watState,ixTopNrg) ! dVol/dT (K-1)
          endif

          ! melt-freeze: compute derivative in energy with respect to mass
          if(mLayerdTheta_dTk(jLayer) > verySmall)then  ! ice is present
            aJac(nrgState,watState) = -dVolTot_dPsi0(iLayer)*LH_fus*iden_water    ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
          else
            aJac(nrgState,watState) = 0._rkind
          endif

          ! - compute lower diagonal elements
          if(iLayer>1)then
            if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer-1),nrgState) = (dt/mLayerDepth(jLayer-1))*( dq_dNrgStateBelow(iLayer-1))   ! K-1
          endif

          ! compute upper-diagonal elements
          if(iLayer<nSoil)then
            if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer+1),nrgState) = (dt/mLayerDepth(jLayer+1))*(-dq_dNrgStateAbove(iLayer))     ! K-1
          endif

        endif   ! (if the energy state for the current layer is within the state subset)

        end do  ! (looping through soil layers)
      endif   ! (if there are state variables for both water and energy in the soil domain)


    ! ***
    ! check
    case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'; return

  end select  ! type of matrix

  ! end association to variables in the data structures
  end associate

end subroutine computJacob


integer(c_int) function computJacob_kinsol(sunvec_y, sunvec_f, sunmat_J, &
                            user_data, sunvec_t1, sunvec_t2 &
                            ) result(ierr) bind(C, name='computJacob_kinsol')
  !======= Inclusions ===========
  USE fsundials_nvector_mod
  USE fsundials_matrix_mod
  USE fnvector_serial_mod
  USE fsunmatrix_dense_mod          
  
  USE kinsol_user_data_type

  implicit none

  ! calling variables
  type(N_Vector)                :: sunvec_y  ! solution N_Vector
  type(N_Vector)                :: sunvec_f  ! rhs N_Vector
  type(SUNMatrix)               :: sunmat_J  ! Jacobian SUNMatrix
  type(c_ptr),value             :: user_data ! user-defined data
  type(N_Vector)                :: sunvec_t1 ! temporary N_Vectors
  type(N_Vector)                :: sunvec_t2

  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer       :: stateVec(:)    ! state vector
  real(c_double), pointer       :: Jac(:,:)       ! Jacobian matrix
  character(len=256)            :: message         ! error message
  type(kinsol_data), pointer    :: kinsol_user_data     ! pointer to user data

  call c_f_pointer(user_data,kinsol_user_data)
  ierr=0

  stateVec => FN_VGetArrayPointer(sunvec_y)
  Jac(1:kinsol_user_data%nState, 1:kinsol_user_data%nstate) => FSUNDenseMatrix_Data(sunmat_J)
  call computJacob(kinsol_user_data%dt,                &
                   kinsol_user_data%nSnow,             &
                   kinsol_user_data%nSoil,             &
                   kinsol_user_data%nLayers,           &
                   kinsol_user_data%computeVegFlux,    &
                   kinsol_user_data%computeBaseflow,   &
                   kinsol_user_data%ixMatrix,          &
                   kinsol_user_data%indx_data,         &
                   kinsol_user_data%prog_data,         &
                   kinsol_user_data%diag_data,         &
                   kinsol_user_data%deriv_data,        &
                   kinsol_user_data%dBaseflow_dMatric, &
                   kinsol_user_data%dMat,              &
                   Jac(:,:),                           &
                   ierr, message)

            
end function computJacob_kinsol

end module computJacob_module
