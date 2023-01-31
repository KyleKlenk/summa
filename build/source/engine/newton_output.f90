! module to output data to netcdf surounding the newton iteration

module newton_output
USE nrtype
USE netcdf

implicit none
private
public::create_newton_file
public::write_linear_system_input
public::write_linear_system_output
! parameter values for the netcdf file
logical(lgt),save,public                    :: newton_file_created=.false. ! flag to indicate if the newton file has been created
integer(i4b),save,public                    :: ncid_newton_file          ! netcdf newton file id
integer(i4b),save,public                    :: ncid_state_dim            ! netcdf state dimension id
integer(i4b),save,public                    :: ncid_iteration_dim        ! netcdf iteration dimension id
integer(i4b),save,public                    :: ncid_nRHS_dim             ! netcdf RHS dimension id       
integer(i4b),save,public                    :: aJac_in_var               ! netcdf aJac_in variable id                   
integer(i4b),save,public                    :: aJac_out_var              ! netcdf aJac_out variable id             
integer(i4b),save,public                    :: ipiv_out_var              ! netcdf ipiv_out variable id            
integer(i4b),save,public                    :: rhs_in_var                ! netcdf rhs_in variable id
integer(i4b),save,public                    :: rhs_out_var               ! netcdf rhs_out variable id
integer(i4b),save,public                    :: iteration_count=1         ! iteration counter
contains
! subroutine to create the netcdf file
subroutine create_newton_file(nState, nRHS)
  implicit none
  integer(i4b),intent(in)     :: nState
  integer(i4b),intent(in)     :: nRHS
  integer(i4b) :: err
  
  ! Define the dimensions
  err = nf90_create(trim("newtonIterationData.nc"),NF90_NETCDF4,ncid_newton_file)
  err = nf90_def_dim(ncid_newton_file, "nState", nState, ncid_state_dim)
  err = nf90_def_dim(ncid_newton_file, "iteration", nf90_unlimited, ncid_iteration_dim)
  err = nf90_def_dim(ncid_newton_file, "nRHS", nRHS, ncid_nRHS_dim)

  ! Define the variables
  err = nf90_def_var(ncid_newton_file, "Jacobian Input for Ax=b", NF90_DOUBLE, (/ncid_iteration_dim, ncid_state_dim, ncid_state_dim/), aJac_in_var)
  err = nf90_def_var(ncid_newton_file, "Jacobian Output for Ax=b", NF90_DOUBLE, (/ncid_iteration_dim, ncid_state_dim, ncid_state_dim/), aJac_out_var)
  err = nf90_def_var(ncid_newton_file, "Pivot Output for Ax=b", NF90_INT, (/ncid_iteration_dim, ncid_state_dim/), ipiv_out_var)
  err = nf90_def_var(ncid_newton_file, "RHS Input for Ax=b", NF90_DOUBLE, (/ncid_iteration_dim, ncid_state_dim, nRHS/), rhs_in_var)
  err = nf90_def_var(ncid_newton_file, "RHS Output for Ax=b", NF90_DOUBLE, (/ncid_iteration_dim, ncid_state_dim, nRHS/), rhs_out_var)

end subroutine create_newton_file

!
subroutine write_linear_system_input(nState, nRHS, jac_input, rhs_input)
  implicit none
  integer(i4b),intent(in)       :: nState
  integer(i4b),intent(in)       :: nRHS
  real(rkind),intent(in)        :: jac_input(:,:)
  real(rkind),intent(in)        :: rhs_input(:,:)

  integer(i4b)                  :: err

  ! Write the data
  err = nf90_put_var(ncid_newton_file, aJac_in_var, jac_input(1:nState, 1:nState),start=(/iteration_count, 1, 1/), count=(/1, nState, nState/))
  err = nf90_put_var(ncid_newton_file, rhs_in_var, rhs_input(1:nState, 1:nRHS),start=(/iteration_count, 1, 1/), count=(/1, nState, nRHS/))

end subroutine write_linear_system_input

subroutine write_linear_system_output(nState, nRHS, jac_output, rhs_output, ipiv_output)
  implicit none
  integer(i4b),intent(in)       :: nState
  integer(i4b),intent(in)       :: nRHS
  real(rkind),intent(in)        :: jac_output(:,:)
  real(rkind),intent(in)        :: rhs_output(:,:)
  integer(i4b),intent(in)       :: ipiv_output(:)


  integer(i4b)                  :: err

  ! Write the data
  err = nf90_put_var(ncid_newton_file, aJac_out_var, jac_output(1:nState, 1:nState),start=(/iteration_count, 1, 1/), count=(/1, nState, nState/))
  err = nf90_put_var(ncid_newton_file, rhs_out_var, rhs_output(1:nState, 1:nRHS),start=(/iteration_count, 1, 1/), count=(/1, nState, nRHS/))
  err = nf90_put_var(ncid_newton_file, ipiv_out_var, ipiv_output(1:nState),start=(/iteration_count, 1/), count=(/1, nState/))

  ! Increment the iteration counter
  iteration_count = iteration_count + 1

end subroutine write_linear_system_output




end module
