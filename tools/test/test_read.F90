!-----------------------------------------------------------------------
!
    program test_cubed_sphere_read
!
!-----------------------------------------------------------------------
!***  This driver is designed to test reading a cubed sphere inc file
!-----------------------------------------------------------------------
!***
!***  Revision history
!***
!     Oct 2022:  M. Potts             - initial code for fv3 read
!
!---------------------------------------------------------------------------------
!
      use mpi, only: MPI_Init, MPI_Finalize
      use module_get_cubed_sphere_inc,  only : read_netcdf
      use CCPP_data,          only: GFS_control
      use fv_arrays_mod,      only: fv_atmos_type
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!
!
      integer              :: rc, grid_id, mype, mpi_comm
      character(128) :: filename
      logical :: tests_passed
      type (fv_atmos_type), allocatable :: Atm(:)
      integer                           :: mygrid

     
      filename = 'atminc.nc4' 
      mpi_comm = 0
      mype = 0
      grid_id = 1
      mygrid = 1
!     call MPI_Init(rc)
      mpi_comm=0
      mype=0
      GFS_control%isc=1
      GFS_control%nx=96
      GFS_control%jsc=1
      GFS_control%ny=96
      call read_netcdf(trim(filename),Atm,mygrid, &
                       .false., &
                       .true., &
        tests_passed=tests_passed,rc=rc)

!     call MPI_Finalize(rc)
      if(tests_passed) then
        call exit(0)
      else
        call exit(1)
      endif
     end program
