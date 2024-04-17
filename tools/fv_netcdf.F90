module module_get_cubed_sphere_inc

  use module_read_netcdf,  only : read_netcdf_inc
  use fv_iau_mod,         only: iau_external_data_type

  implicit none


  contains

  subroutine get_cubed_sphere_inc(filename, increment_data, &
                          testing, tests_passed,rc)
!
    character(*), intent(in)           :: filename
    type( IAU_external_data_type )     :: increment_data
    logical, intent(in)                :: testing
    logical, optional,intent(out)      :: tests_passed
    integer, optional,intent(out)      :: rc
!
      
  end subroutine get_cubed_sphere_inc
  
!----------------------------------------------------------------------------------------
end module module_get_cubed_sphere_inc
