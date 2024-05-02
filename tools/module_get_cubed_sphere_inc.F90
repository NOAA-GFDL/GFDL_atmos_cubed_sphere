module module_get_cubed_sphere_inc

  use netcdf
  use esmf
  use mpp_mod,            only: mpp_pe, mpp_get_current_pelist
  use tracer_manager_mod, only: get_tracer_index, get_number_tracers, get_tracer_names
  use field_manager_mod,  only: MODEL_ATMOS
  use CCPP_data,          only: GFS_control
  use time_manager_mod,   only: time_type, get_time, set_time
  use fv_arrays_mod,      only: fv_atmos_type
  use fms2_io_mod,        only: open_file, close_file, read_data, variable_exists, FmsNetcdfFile_t
  use mpi

  implicit none
  type iau_internal_data_type
    real, allocatable :: ua_inc(:,:,:)
    real, allocatable :: va_inc(:,:,:)
    real, allocatable :: temp_inc(:,:,:)
    real, allocatable :: delp_inc(:,:,:)
    real, allocatable :: delz_inc(:,:,:)
    real, allocatable :: tracer_inc(:,:,:,:)
  end type iau_internal_data_type

  public read_netcdf_inc, iau_internal_data_type

  contains

!----------------------------------------------------------------------------------------

    subroutine read_netcdf_inc(filename, increment_data, Atm)
      character(*), intent(in)                    :: filename
      type(iau_internal_data_type), intent(inout) :: increment_data
      type(fv_atmos_type), intent(in)             :: Atm

      type(FmsNetcdfFile_t) :: fileobj
      integer               :: isc, iec, jsc, jec, npz, itracer, ntracers, itile
      integer, dimension(4) :: corner, edge_lengths
      character(len=64)     :: tracer_name
      
      ! Get various dimensions
      call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
      isc = Atm%bd%isc
      iec = Atm%bd%iec
      jsc = Atm%bd%jsc
      jec = Atm%bd%jec
      npz = Atm%npz
      itile = Atm%tile_of_mosaic

      ! Open increment file
      if ( open_file(fileobj, trim(filename), 'read', pelist=Atm%pelist) ) then
         
         !
         corner       = (/ isc,       jsc,       1,   itile /)
         edge_lengths = (/ iec-isc+1, jec-jsc+1, npz, 1 /)

         ! Read increments
         call read_data(fileobj, 'u_inc',    increment_data%ua_inc,   corner=corner, edge_lengths=edge_lengths)
         call read_data(fileobj, 'v_inc',    increment_data%va_inc,   corner=corner, edge_lengths=edge_lengths)
         call read_data(fileobj, 'T_inc',    increment_data%temp_inc, corner=corner, edge_lengths=edge_lengths)
         call read_data(fileobj, 'delp_inc', increment_data%delp_inc, corner=corner, edge_lengths=edge_lengths)
         if ( .not. Atm%flagstruct%hydrostatic ) then
            call read_data(fileobj, 'delz_inc', increment_data%delz_inc, corner=corner, edge_lengths=edge_lengths)
         endif

         ! Read tracer increments
         do itracer = 1,ntracers
            call get_tracer_names(MODEL_ATMOS, itracer, tracer_name)
            if ( variable_exists(fileobj, trim(tracer_name)//'_inc') ) then
               call read_data(fileobj, trim(tracer_name)//'_inc', increment_data%tracer_inc(:,:,:,itracer), &
                    corner=corner, edge_lengths=edge_lengths)
            end if
         end do

         ! Close increment file
         call close_file(fileobj)
         
      end if
    
  end subroutine read_netcdf_inc
    
!----------------------------------------------------------------------------------------
end module module_get_cubed_sphere_inc
