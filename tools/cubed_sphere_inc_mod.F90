module cubed_sphere_inc_mod

  use tracer_manager_mod, only: get_tracer_index, get_number_tracers, get_tracer_names
  use field_manager_mod,  only: MODEL_ATMOS
  use fv_arrays_mod,      only: fv_atmos_type
  use fms2_io_mod,        only: open_file, close_file, read_data, variable_exists, FmsNetcdfDomainFile_t
  
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

    subroutine read_cubed_sphere_inc(fname, increment_data, Atm)
      character(*), intent(in)                 :: fname
      type(increment_data_type), intent(inout) :: increment_data
      type(fv_atmos_type), intent(in)          :: Atm

      type(FmsNetcdfDomainFile_t) :: fileobj
      integer                     :: itracer, ntracers
      character(len=64)           :: tracer_name
      
      ! Get various dimensions
      call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)

      ! Open file
      if ( open_file(fileobj, trim(fname), 'read', Atm%domain) ) then
         ! Read increments
         call read_data(fileobj, 'u_inc', increment_data%ua_inc)
         call read_data(fileobj, 'v_inc', increment_data%va_inc)
         call read_data(fileobj, 'T_inc', increment_data%temp_inc)
         call read_data(fileobj, 'delp_inc', increment_data%delp_inc)
         if ( .not. Atm%flagstruct%hydrostatic ) then
            call read_data(fileobj, 'delz_inc', increment_data%delz_inc)
         end if

         ! Read tracer increments
         do itracer = 1,ntracers
            call get_tracer_names(MODEL_ATMOS, itracer, tracer_name)
            if ( variable_exists(fileobj, trim(tracer_name)//'_inc') ) then
               call read_data(fileobj, trim(tracer_name)//'_inc', increment_data%tracer_inc(:,:,:,itracer))
            end if
         end do
         
         call close_file(fileobj)
      end if
    
  end subroutine read_cubed_sphere_inc

!----------------------------------------------------------------------------------------
end module cubed_sphere_inc_mod
