module cubed_sphere_inc_mod

  use tracer_manager_mod, only: get_tracer_index, get_number_tracers, get_tracer_names
  use field_manager_mod,  only: MODEL_ATMOS
  use fv_arrays_mod,      only: fv_atmos_type
  use sim_nc_mod,         only: open_ncfile, get_var5_r4, check_var_exists, close_ncfile

  implicit none
  type increment_data_type
    real, allocatable :: ua_inc(:,:,:)
    real, allocatable :: va_inc(:,:,:)
    real, allocatable :: temp_inc(:,:,:)
    real, allocatable :: delp_inc(:,:,:)
    real, allocatable :: delz_inc(:,:,:)
    real, allocatable :: tracer_inc(:,:,:,:)
  end type increment_data_type

  public read_cubed_sphere_inc, increment_data_type

  contains

!----------------------------------------------------------------------------------------

    subroutine read_cubed_sphere_inc(filename, increment_data, Atm)
      character(*), intent(in)                 :: filename
      type(increment_data_type), intent(inout) :: increment_data
      type(fv_atmos_type), intent(in)          :: Atm

      integer           :: isc, iec, jsc, jec, npz, itracer, ntracers, itile
      character(len=64) :: tracer_name
      integer           :: ncid, ncerr
      
      ! Get various dimensions
      call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
      isc = Atm%bd%isc
      iec = Atm%bd%iec
      jsc = Atm%bd%jsc
      jec = Atm%bd%jec
      npz = Atm%npz
      itile = Atm%tile_of_mosaic

      ! Open increment file
      call open_ncfile( trim(filename), ncid )

      ! Read increments
      call get_var5_r4( ncid, 'u_inc', isc,iec, jsc,jec, 1,npz, itile, 1, increment_data%ua_inc )
      call get_var5_r4( ncid, 'v_inc', isc,iec, jsc,jec, 1,npz, itile, 1, increment_data%va_inc )
      call get_var5_r4( ncid, 'temp_inc', isc,iec, jsc,jec, 1,npz, itile, 1, increment_data%temp_inc )
      call get_var5_r4( ncid, 'delp_inc', isc,iec, jsc,jec, 1,npz, itile, 1, increment_data%delp_inc )
      if ( .not. Atm%flagstruct%hydrostatic ) then
         call get_var5_r4( ncid, 'delz_inc', isc,iec, jsc,jec, 1,npz, itile, 1, increment_data%delz_inc )
      end if

      ! Read tracer increments
      do itracer = 1,ntracers
         call get_tracer_names(MODEL_ATMOS, itracer, tracer_name)
         call check_var_exists(ncid, trim(tracer_name)//'_inc', ncerr)
         if ( ncerr == 0 ) then
            call get_var5_r4( ncid, trim(tracer_name)//'_inc', isc,iec, jsc,jec, 1,npz, itile, 1, increment_data%tracer_inc(:,:,:,itracer) )
         end if
      end do
         
      ! Close increment file
      call close_ncfile(ncid)
    
  end subroutine read_cubed_sphere_inc

!----------------------------------------------------------------------------------------
end module cubed_sphere_inc_mod
