!-----------------------------------------------------------------------
!
    program test_iau
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
      use atmosphere_mod,     only: Atm, mygrid

      use module_get_cubed_sphere_inc,  only : calc_iau_tendency, iau_internal_data_type
      use CCPP_data,          only: GFS_control
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
      character(len=128), dimension(3) :: filenames
      type(iau_internal_data_type)     :: tendency ! number of blocks
      type(iau_internal_data_type)     :: increment ! number of blocks
!     type(GFS_init_type)  :: Init_parm
      logical :: tests_passed, use_parallel_netcdf, testing
      integer :: itest, jtest, ktest, i
      real :: u, dt_atmos, incsum, incnum

      integer :: n, clock
      integer, dimension(:), allocatable :: seed
          
      call random_seed(size = n)
      allocate(seed(n))
          
      call system_clock(count=clock)
          
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)
          
      deallocate(seed) 
      filenames = (/'temp0.nc4', 'temp1.nc4', 'temp2.nc4'/)
      use_parallel_netcdf=.false.
      testing=.true.
      mpi_comm = 0
      mype = 0
      grid_id = 1
      mpi_comm=0
      mype=0
      dt_atmos = 250.
      GFS_control%isc=1
      GFS_control%nx=96
      GFS_control%jsc=1
      GFS_control%ny=96
      call calc_iau_tendency(filenames, Atm, mygrid, &
                       tendency, increment=increment, &
                       testing=testing, rc=rc)

      incsum = 0.0
      incnum = 0
      ! The increment files were generated with a perterbation of 0.01 up and down from temp1.nc4
      ! to test that we are getting values read in correctly, we will do a bit of algebra to 
      ! make sure the differences are all 0.01 (or average very closely out to that)
      do i = 1,10
        call random_number(u)
        itest = floor(96*u) + 1
        call random_number(u)
        jtest = floor(96*u) + 1
        call random_number(u)
        ktest = floor(127*u) + 1
    
        if(increment%ua_inc(itest,jtest,ktest) /= 0.0) then
           write(6,*) tendency%ua_inc(itest,jtest,ktest)*dt_atmos/increment%ua_inc(itest,jtest,ktest)
           incsum = incsum + tendency%ua_inc(itest,jtest,ktest)*dt_atmos/increment%ua_inc(itest,jtest,ktest)
           incnum = incnum + 1
        endif
        if(increment%va_inc(itest,jtest,ktest) /= 0.0) then
           incsum = incsum + tendency%va_inc(itest,jtest,ktest)*dt_atmos/increment%va_inc(itest,jtest,ktest)
           incnum = incnum + 1
        endif
        if(increment%temp_inc(itest,jtest,ktest) /= 0.0) then
           incsum = incsum + tendency%temp_inc(itest,jtest,ktest)*dt_atmos/increment%temp_inc(itest,jtest,ktest)
           incnum = incnum + 1
        endif
        if(increment%delp_inc(itest,jtest,ktest) /=  0.0) then
           incsum = incsum + tendency%delp_inc(itest,jtest,ktest)*dt_atmos/increment%delp_inc(itest,jtest,ktest)
           incnum = incnum + 1
        endif
        if(increment%tracer_inc(itest,jtest,ktest,1) > 0.0) then
           incsum = incsum + tendency%tracer_inc(itest,jtest,ktest,1)*dt_atmos/increment%tracer_inc(itest,jtest,ktest,1)
           incnum = incnum + 1
        endif
        if(increment%tracer_inc(itest,jtest,ktest,2) > 0.0) then
           incsum = incsum + tendency%tracer_inc(itest,jtest,ktest,2)*dt_atmos/increment%tracer_inc(itest,jtest,ktest,2)
           incnum = incnum + 1
        endif
        if(increment%tracer_inc(itest,jtest,ktest,3) > 0.0) then
           incsum = incsum + tendency%tracer_inc(itest,jtest,ktest,3)*dt_atmos/increment%tracer_inc(itest,jtest,ktest,3)
           incnum = incnum + 1
        endif
        if(increment%tracer_inc(itest,jtest,ktest,4) > 0.0) then
           incsum = incsum + tendency%tracer_inc(itest,jtest,ktest,4)*dt_atmos/increment%tracer_inc(itest,jtest,ktest,4)
           incnum = incnum + 1
        endif
      enddo
      write(6,*) 'incsum, incnum and diffs is ',incsum,incnum,incnum * 0.01 
      if( abs( incsum - incnum * 0.01 ) < 0.000001) then
         tests_passed = .true.
      else
         tests_passed = .false.
      endif
      if(tests_passed) then
        call exit(0)
      else
        call exit(1)
      endif
     end program
