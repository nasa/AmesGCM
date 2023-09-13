module aerosol_util_mod

use fms_mod,                only: error_mesg, FATAL,       &
                                 open_namelist_file, check_nml_error, &
                                 mpp_pe, mpp_root_pe, close_file,     &
                                 write_version_number, stdlog,        &
                                 uppercase, read_data, write_data, field_size

use fms2_io_mod,            only:  file_exists

implicit none
private

public :: init_aerosol_flags
public :: dust_map_scale, dust_map_scale_bin

logical, public ::  do_moment_dust  = .true.             ! do moment dust lifting
logical, public ::  do_moment_water = .false.            ! do moment water microphysics
logical, public ::  do_moment_sedim = .true.             ! do moment sedimentation
logical, public ::  do_15band = .false.
real, public    ::  Reff_backgd = 2.0e-6                 !  effective radius for lifting in the case of background scenario
real, public    ::  Reff_stress = 2.0e-6                 !  effective radius for stress lifting
real, public    ::  Reff_dd = 2.0e-6                     !  effective radius for dust devils lifting
real :: dust_map_scale = 3.67
real :: dust_map_scale_bin = 3.67
real, public    ::  Reff_fixed = 1.5e-6                  !  effective radius for fixed background dust

namelist /aerosol_util_nml/  do_moment_dust, do_moment_water, dust_map_scale, dust_map_scale_bin, &
                             Reff_backgd, Reff_stress, Reff_dd, do_moment_sedim, do_15band, &
                             Reff_fixed


contains


subroutine init_aerosol_flags

integer  unit, io, ierr


!     ----- read namelist -----

if (file_exists('input.nml')) then
    unit = open_namelist_file ( )
    ierr=1; do while (ierr /= 0)
        read  (unit, nml=aerosol_util_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'aerosol_util_nml')
    enddo
10     call close_file (unit)
endif

if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=aerosol_util_nml)

if (do_moment_water .and. .not. do_moment_dust) call error_mesg ('aerosol_util','moment water microphysics turned on without moment dust', FATAL)

end subroutine init_aerosol_flags


end module aerosol_util_mod
