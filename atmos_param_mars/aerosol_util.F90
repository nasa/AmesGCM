module aerosol_util_mod


use           fms_mod, only: error_mesg, FATAL,       &
                             check_nml_error, &
                             mpp_pe, mpp_root_pe, &
                             write_version_number, stdlog,        &
                             uppercase

use       fms2_io_mod, only:  file_exists, FmsNetcdfFile_t, FmsNetcdfDomainFile_t, &
                                   register_restart_field, register_axis, unlimited, &
                                   open_file, read_restart, write_restart, close_file, &
                                   register_field, read_data, write_data, register_variable_attribute, &
                                   get_global_io_domain_indices, get_variable_size, variable_exists

use   mpp_mod, only: input_nml_file
implicit none
private

public :: init_aerosol_flags
public :: dust_map_scale, dust_map_scale_bin

logical, public ::  do_moment_dust  = .true.             ! do moment dust lifting
logical, public ::  do_moment_water = .false.            ! do moment water microphysics
logical, public ::  do_moment_sedim = .true.             ! do moment sedimentation
logical, public ::  do_bulk_water = .false.              ! do bulk water cloud physics
logical, public ::  do_15band = .false.
real, public    ::  Reff_backgd = 2.0e-6                 !  effective radius for lifting in the case of background scenario
real, public    ::  Reff_stress = 2.0e-6                 !  effective radius for stress lifting
real, public    ::  Reff_dd = 2.0e-6                     !  effective radius for dust devils lifting
real :: dust_map_scale = 3.67
real :: dust_map_scale_bin = 3.67
real, public    ::  Reff_fixed = 1.5e-6                  !  effective radius for fixed background dust

namelist /aerosol_util_nml/  do_moment_dust, do_moment_water, dust_map_scale, dust_map_scale_bin, &
                             do_bulk_water, &
                             Reff_backgd, Reff_stress, Reff_dd, do_moment_sedim, do_15band, &
                             Reff_fixed


contains


subroutine init_aerosol_flags

integer  unit, io, ierr


!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
read (input_nml_file, nml=aerosol_util_nml, iostat=io)
ierr = check_nml_error(io,'aerosol_util_nml')

if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=aerosol_util_nml)

if (do_moment_water .and. .not. do_moment_dust) call error_mesg ('aerosol_util','moment water microphysics turned on without moment dust', FATAL)

if (do_moment_water .and. do_bulk_water) call error_mesg ('aerosol_util','moment water microphysics and bulk water physics cannot both be turned on', FATAL)

if (do_moment_dust .and. do_bulk_water) call error_mesg ('aerosol_util','moment dust microphysics and bulk water physics cannot both be turned on', FATAL)

end subroutine init_aerosol_flags


end module aerosol_util_mod
