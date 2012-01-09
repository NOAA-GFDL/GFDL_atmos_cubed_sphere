!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fv_mp_mod --- SPMD parallel decompostion/communication module
      module fv_mp_mod

#if defined(SPMD)
! !USES:
      use fms_mod,         only : fms_init, fms_end
      use mpp_mod,         only : FATAL, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
      use mpp_mod,         only : mpp_pe, mpp_npes, mpp_node, mpp_root_pe, mpp_error, mpp_set_warn_level
      use mpp_mod,         only : mpp_declare_pelist, mpp_set_current_pelist, mpp_sync
      use mpp_mod,         only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
      use mpp_mod,         only : mpp_chksum, stdout, stderr, mpp_broadcast
      use mpp_mod,         only : mpp_send, mpp_recv, mpp_sync_self, EVENT_RECV, mpp_gather
      use mpp_domains_mod, only : GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, FOLD_NORTH_EDGE, CGRID_NE
      use mpp_domains_mod, only : MPP_DOMAIN_TIME, CYCLIC_GLOBAL_DOMAIN, NUPDATE,EUPDATE, XUPDATE, YUPDATE, SCALAR_PAIR
      use mpp_domains_mod, only : domain1D, domain2D, DomainCommunicator2D
      use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
      use mpp_domains_mod, only : mpp_global_field, mpp_global_sum, mpp_global_max, mpp_global_min
      use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
      use mpp_domains_mod, only : mpp_check_field, mpp_define_layout 
      use mpp_domains_mod, only : mpp_get_neighbor_pe, mpp_define_mosaic, mpp_define_io_domain
      use mpp_domains_mod, only : NORTH, NORTH_EAST, EAST, SOUTH_EAST
      use mpp_domains_mod, only : SOUTH, SOUTH_WEST, WEST, NORTH_WEST
      use mpp_parameter_mod, only : WUPDATE, EUPDATE, SUPDATE, NUPDATE, XUPDATE, YUPDATE

      implicit none
      private

#include "mpif.h"
      integer, parameter :: XDir=1
      integer, parameter :: YDir=2
      integer :: npes, npes_x, npes_y, gid, masterproc, commglobal, ierror
      integer :: io_domain_layout(2) = (/1,1/)

      type(domain2D), target, save :: domain
      type(domain2D), target, save :: domain_for_coupler ! domain used in coupled model with halo = 1.
      integer :: num_contact, ntiles, npes_per_tile, tile
      integer, allocatable, dimension(:)       :: npes_tile, tile1, tile2
      integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
      integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
      integer, allocatable, dimension(:,:)     :: layout2D, global_indices
      integer, parameter:: ng    = 3     ! Number of ghost zones required
      integer :: is, ie, js, je, isd, ied, jsd, jed
      integer :: numthreads
      logical:: square_domain = .false.

      
      public mp_start, mp_barrier, mp_stop, npes, npes_x, npes_y, ng, gid, masterproc
      public io_domain_layout, square_domain
      public domain, tile, domain_for_coupler
      public is, ie, js, je, isd, ied, jsd, jed
      public domain_decomp, mp_bcst, mp_reduce_max, mp_reduce_sum, mp_gather
      public mp_reduce_min
      public fill_corners, mp_corner_comm, XDir, YDir
      public numthreads

      INTERFACE fill_corners
        MODULE PROCEDURE fill_corners_2d
        MODULE PROCEDURE fill_corners_xy_2d
        MODULE PROCEDURE fill_corners_xy_3d
      END INTERFACE

      INTERFACE mp_bcst
        MODULE PROCEDURE mp_bcst_i4
        MODULE PROCEDURE mp_bcst_r8
        MODULE PROCEDURE mp_bcst_3d_r8
        MODULE PROCEDURE mp_bcst_4d_r8
      END INTERFACE

      INTERFACE mp_reduce_max
        MODULE PROCEDURE mp_reduce_max_r8_1d
        MODULE PROCEDURE mp_reduce_max_r8
        MODULE PROCEDURE mp_reduce_max_i4
      END INTERFACE

      INTERFACE mp_reduce_sum
        MODULE PROCEDURE mp_reduce_sum_r8
        MODULE PROCEDURE mp_reduce_sum_r8_1d
      END INTERFACE

      INTERFACE mp_gather
        MODULE PROCEDURE mp_gather_4d_r4
        MODULE PROCEDURE mp_gather_3d_r4
        MODULE PROCEDURE mp_gather_3d_r8
      END INTERFACE

!---- version number -----
      character(len=128) :: version = '$Id: fv_mp_mod.F90,v 19.0 2012/01/06 19:59:11 fms Exp $'
      character(len=128) :: tagname = '$Name: siena $'

contains

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_start :: Start SPMD processes
!
        subroutine mp_start(commID)
         integer, intent(in), optional :: commID
         integer :: unit
!$       integer :: omp_get_num_threads

         gid = mpp_pe()
         npes = mpp_npes()

         masterproc = mpp_root_pe()
         commglobal = MPI_COMM_WORLD
         if( PRESENT(commID) )then
             commglobal = commID
         end if

         numthreads = 1
!$OMP PARALLEL
!$OMP MASTER
!$       numthreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
         if ( mpp_pe()==mpp_root_pe() ) then
           unit = stdout()
           write(unit,*) 'Starting PEs : ', npes
           write(unit,*) 'Starting Threads : ', numthreads
         endif

         call MPI_BARRIER(commglobal, ierror)
      end subroutine mp_start
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_barrier :: Wait for all SPMD processes
!
      subroutine mp_barrier()
        
         call MPI_BARRIER(commglobal, ierror)
      
      end subroutine mp_barrier
!       
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_stop :: Stop SPMD processes
!
      subroutine mp_stop()

         call MPI_BARRIER(commglobal, ierror)
         if (gid==masterproc) print*, 'Stopping PEs : ', npes
         call fms_end()
        ! call MPI_FINALIZE (ierror)

      end subroutine mp_stop
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     domain_decomp :: Setup domain decomp
!
      subroutine domain_decomp(npx,npy,nregions,ng,grid_type)

         integer, intent(IN)  :: npx,npy,nregions,ng,grid_type

         integer :: layout(2)
         integer, allocatable :: pe_start(:), pe_end(:)

         character(len=80) :: evalue
         integer :: ios,nx,ny,n,num_alloc
         character(len=32) :: type = "unknown"
         logical :: is_symmetry 
         logical :: debug=.false.

         nx = npx-1
         ny = npy-1

         call mpp_domains_init(MPP_DOMAIN_TIME)

       ! call mpp_domains_set_stack_size(10000)
       ! call mpp_domains_set_stack_size(900000)
         call mpp_domains_set_stack_size(1500000)

         select case(nregions)
         case ( 1 )  ! Lat-Lon "cyclic"

            select case (grid_type)
            case (3)   ! Lat-Lon "cyclic"
               type="Lat-Lon: cyclic"
               ntiles = 4
               num_contact = 8
               if( mod(npes,ntiles) .NE. 0 ) then
                  call mpp_error(NOTE,'TEST_MPP_DOMAINS: for Cyclic mosaic, npes should be multiple of ntiles. ' // &
                                       'No test is done for Cyclic mosaic. ' )
                  return
               end if
               npes_per_tile = npes/ntiles
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
               layout = (/1,npes_per_tile/) ! force decomp only in Lat-Direction
            case (4)   ! Cartesian, double periodic
               type="Cartesian: double periodic"
               ntiles = 1
               num_contact = 2
               npes_per_tile = npes/ntiles
               if(npes_x*npes_y == npes_per_tile) then
                  layout = (/npes_x,npes_y/)
               else
                  call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
               endif
            case (5)   ! latlon patch
               type="Lat-Lon: patch"
               ntiles = 1
               num_contact = 0
               npes_per_tile = npes/ntiles
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
            case (6)   ! latlon strip
               type="Lat-Lon: strip"
               ntiles = 1
               num_contact = 1
               npes_per_tile = npes/ntiles
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
            case (7)   ! Cartesian, channel
               type="Cartesian: channel"
               ntiles = 1
               num_contact = 1
               npes_per_tile = npes/ntiles
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
            end select

         case ( 6 )  ! Cubed-Sphere
            type="Cubic: cubed-sphere"
            ntiles = 6
            num_contact = 12
            !--- cubic grid always have six tiles, so npes should be multiple of 6
            if( mod(npes,ntiles) .NE. 0 .OR. npx-1 .NE. npy-1) then
               call mpp_error(NOTE,'domain_decomp: for Cubic_grid mosaic, npes should be multiple of ntiles(6) ' // &
                                   'and npx-1 should equal npy-1, domain_decomp is NOT done for Cubic-grid mosaic. ' )
               return
            end if
            npes_per_tile = npes/ntiles
            call  mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )

            if ( npes_x == 0 ) then 
               npes_x = layout(1)
            endif
            if ( npes_y == 0 ) then
               npes_y = layout(2)
            endif

            if ( npes_x==npes_y .and. (npx-1)==((npx-1)/npes_x)*npes_x )  square_domain = .true.

            if ( (npx/npes_x < ng) .or. (npy/npes_y < ng) ) then
               write(*,310) npes_x, npes_y, npx/npes_x, npy/npes_y
 310           format('Invalid layout, NPES_X:',i4.4,'NPES_Y:',i4.4,'ncells_X:',i4.4,'ncells_Y:',i4.4)
               call mp_stop
               call exit(1)
            endif
           
            layout = (/npes_x,npes_y/)
         case default
            call mpp_error(FATAL, 'domain_decomp: no such test: '//type)
         end select

         allocate(layout2D(2,ntiles), global_indices(4,ntiles), npes_tile(ntiles) )
         allocate(pe_start(ntiles),pe_end(ntiles))
         npes_tile = npes_per_tile
         do n = 1, ntiles
            global_indices(:,n) = (/1,npx-1,1,npy-1/)
            layout2D(:,n)         = layout
            pe_start(n) = mpp_root_pe() + (n-1)*layout(1)*layout(2)
            pe_end(n)   = pe_start(n) + layout(1)*layout(2) -1
         end do
         num_alloc=max(1,num_contact)
         allocate(tile1(num_alloc), tile2(num_alloc) )
         allocate(istart1(num_alloc), iend1(num_alloc), jstart1(num_alloc), jend1(num_alloc) )
         allocate(istart2(num_alloc), iend2(num_alloc), jstart2(num_alloc), jend2(num_alloc) )
 
         is_symmetry = .true.
         select case(nregions)
         case ( 1 )

            select case (grid_type)
            case (3)   ! Lat-Lon "cyclic"
               !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
               tile1(1) = 1; tile2(1) = 2
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
               !--- Contact line 2, between tile 1 (SOUTH) and tile 3 (NORTH)  --- cyclic
               tile1(2) = 1; tile2(2) = 3
               istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
               istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
               !--- Contact line 3, between tile 1 (WEST) and tile 2 (EAST) --- cyclic
               tile1(3) = 1; tile2(3) = 2
               istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
               istart2(3) = nx; iend2(3) = nx; jstart2(3) = 1;  jend2(3) = ny
               !--- Contact line 4, between tile 1 (NORTH) and tile 3 (SOUTH)
               tile1(4) = 1; tile2(4) = 3
               istart1(4) = 1;  iend1(4) = nx; jstart1(4) = ny;  jend1(4) = ny
               istart2(4) = 1;  iend2(4) = nx; jstart2(4) = 1;   jend2(4) = 1
               !--- Contact line 5, between tile 2 (SOUTH) and tile 4 (NORTH) --- cyclic
               tile1(5) = 2; tile2(5) = 4
               istart1(5) = 1;  iend1(5) = nx; jstart1(5) = 1;  jend1(5) = 1
               istart2(5) = 1;  iend2(5) = nx; jstart2(5) = ny; jend2(5) = ny
               !--- Contact line 6, between tile 2 (NORTH) and tile 4 (SOUTH)
               tile1(6) = 2; tile2(6) = 4
               istart1(6) = 1;  iend1(6) = nx; jstart1(6) = ny;  jend1(6) = ny
               istart2(6) = 1;  iend2(6) = nx; jstart2(6) = 1;   jend2(6) = 1
               !--- Contact line 7, between tile 3 (EAST) and tile 4 (WEST)
               tile1(7) = 3; tile2(7) = 4
               istart1(7) = nx; iend1(7) = nx; jstart1(7) = 1;  jend1(7) = ny
               istart2(7) = 1;  iend2(7) = 1;  jstart2(7) = 1;  jend2(7) = ny
               !--- Contact line 8, between tile 3 (WEST) and tile 4 (EAST) --- cyclic
               tile1(8) = 3; tile2(8) = 4
               istart1(8) = 1;  iend1(8) = 1;  jstart1(8) = 1;  jend1(8) = ny
               istart2(8) = nx; iend2(8) = nx; jstart2(8) = 1;  jend2(8) = ny
               is_symmetry = .false.
            case (4)   ! Cartesian, double periodic
               !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
               tile1(1) = 1; tile2(1) = 1
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
               !--- Contact line 2, between tile 1 (SOUTH) and tile 1 (NORTH)  --- cyclic
               tile1(2) = 1; tile2(2) = 1
               istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
               istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
            case (5)   ! latlon patch

            case (6)   !latlon strip
               !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
               tile1(1) = 1; tile2(1) = 1
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
            case (7)   ! Cartesian, channel
               !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
               tile1(1) = 1; tile2(1) = 1
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
            end select

         case ( 6 )  ! Cubed-Sphere
            !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
            tile1(1) = 1; tile2(1) = 2
            istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
            istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
            !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
            tile1(2) = 1; tile2(2) = 3
            istart1(2) = 1;  iend1(2) = nx; jstart1(2) = ny; jend1(2) = ny
            istart2(2) = 1;  iend2(2) = 1;  jstart2(2) = ny; jend2(2) = 1
            !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
            tile1(3) = 1; tile2(3) = 5
            istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
            istart2(3) = nx; iend2(3) = 1;  jstart2(3) = ny; jend2(3) = ny
            !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
            tile1(4) = 1; tile2(4) = 6
            istart1(4) = 1;  iend1(4) = nx; jstart1(4) = 1;  jend1(4) = 1
            istart2(4) = 1;  iend2(4) = nx; jstart2(4) = ny; jend2(4) = ny
            !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
            tile1(5) = 2; tile2(5) = 3
            istart1(5) = 1;  iend1(5) = nx; jstart1(5) = ny; jend1(5) = ny
            istart2(5) = 1;  iend2(5) = nx; jstart2(5) = 1;  jend2(5) = 1
            !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
            tile1(6) = 2; tile2(6) = 4
            istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
            istart2(6) = nx; iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = 1
            !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
            tile1(7) = 2; tile2(7) = 6
            istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;  jend1(7) = 1
            istart2(7) = nx; iend2(7) = nx; jstart2(7) = ny; jend2(7) = 1
            !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
            tile1(8) = 3; tile2(8) = 4
            istart1(8) = nx; iend1(8) = nx; jstart1(8) = 1;  jend1(8) = ny
            istart2(8) = 1;  iend2(8) = 1;  jstart2(8) = 1;  jend2(8) = ny
            !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
            tile1(9) = 3; tile2(9) = 5
            istart1(9) = 1;  iend1(9) = nx; jstart1(9) = ny; jend1(9) = ny
            istart2(9) = 1;  iend2(9) = 1;  jstart2(9) = ny; jend2(9) = 1
            !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
            tile1(10) = 4; tile2(10) = 5
            istart1(10) = 1;  iend1(10) = nx; jstart1(10) = ny; jend1(10) = ny
            istart2(10) = 1;  iend2(10) = nx; jstart2(10) = 1;  jend2(10) = 1
            !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
            tile1(11) = 4; tile2(11) = 6
            istart1(11) = nx; iend1(11) = nx; jstart1(11) = 1;  jend1(11) = ny
            istart2(11) = nx; iend2(11) = 1;  jstart2(11) = 1;  jend2(11) = 1
            !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
            tile1(12) = 5; tile2(12) = 6
            istart1(12) = nx; iend1(12) = nx; jstart1(12) = 1;  jend1(12) = ny
            istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;  jend2(12) = ny
         end select

       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              pe_start=pe_start, pe_end=pe_end, symmetry=is_symmetry,              &
                              shalo = ng, nhalo = ng, whalo = ng, ehalo = ng, name = type)
       call mpp_define_mosaic(global_indices, layout2D, domain_for_coupler, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,                  &
                              pe_start=pe_start, pe_end=pe_end, symmetry=is_symmetry,                          &
                              shalo = 1, nhalo = 1, whalo = 1, ehalo = 1, name = type)

       call mpp_define_io_domain(domain, io_domain_layout)
       call mpp_define_io_domain(domain_for_coupler, io_domain_layout)

       deallocate(pe_start,pe_end)

        !--- find the tile number
         tile = (mpp_pe()-mpp_root_pe())/npes_per_tile+1
         call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
         call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

         if (debug .and. nregions==1) then
            tile=1
            write(*,200) tile, is, ie, js, je
         !   call mp_stop
         !   stop
         endif
 200     format(i4.4, ' ', i4.4, ' ', i4.4, ' ', i4.4, ' ', i4.4, ' ')

      end subroutine domain_decomp
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !!
!     
      subroutine fill_corners_2d(q, npx, npy, FILL, AGRID, BGRID)
         real, DIMENSION(isd:,jsd:), intent(INOUT):: q
         integer, intent(IN):: npx,npy
         integer, intent(IN):: FILL  ! X-Dir or Y-Dir 
         logical, OPTIONAL, intent(IN) :: AGRID, BGRID 
         integer :: i,j

         if (present(BGRID)) then
            if (BGRID) then
              select case (FILL)
              case (XDir)
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-i  ,1-j  ) = q(1-j  ,i+1    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-i  ,npy+j) = q(1-j  ,npy-i  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+i,1-j  ) = q(npx+j,i+1    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+i,npy+j) = q(npx+j,npy-i  )  !NE Corner
                    enddo
                 enddo
              case (YDir)
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-j  ,1-i  ) = q(i+1  ,1-j    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-j  ,npy+i) = q(i+1  ,npy+j  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+j,1-i  ) = q(npx-i,1-j    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+j,npy+i) = q(npx-i,npy+j  )  !NE Corner
                    enddo
                 enddo
              case default
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-i  ,1-j  ) = q(1-j  ,i+1    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-i  ,npy+j) = q(1-j  ,npy-i  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+i,1-j  ) = q(npx+j,i+1    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+i,npy+j) = q(npx+j,npy-i  )  !NE Corner
                    enddo
                 enddo
              end select
            endif
          elseif (present(AGRID)) then
            if (AGRID) then
              select case (FILL)
              case (XDir)
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) q(1-i    ,1-j    ) = q(1-j    ,i        )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-i    ,npy-1+j) = q(1-j    ,npy-1-i+1)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+i,1-j    ) = q(npx-1+j,i        )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+i,npy-1+j) = q(npx-1+j,npy-1-i+1)  !NE Corner
                    enddo
                 enddo
              case (YDir)
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) q(1-j    ,1-i    ) = q(i        ,1-j    )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-j    ,npy-1+i) = q(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+j,1-i    ) = q(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+j,npy-1+i) = q(npx-1-i+1,npy-1+j)  !NE Corner
                    enddo
                 enddo
              case default
                 do j=1,ng
                    do i=1,ng        
                       if ((is==    1) .and. (js==    1)) q(1-j    ,1-i    ) = q(i        ,1-j    )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-j    ,npy-1+i) = q(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+j,1-i    ) = q(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+j,npy-1+i) = q(npx-1-i+1,npy-1+j)  !NE Corner
                   enddo
                 enddo          
              end select
            endif
          endif

      end subroutine fill_corners_2d
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !!
!     fill_corners_xy_2d
      subroutine fill_corners_xy_2d(x, y, npx, npy, DGRID, AGRID, CGRID, VECTOR)
         real, DIMENSION(isd:,jsd:), intent(INOUT):: x !(isd:ied  ,jsd:jed+1)
         real, DIMENSION(isd:,jsd:), intent(INOUT):: y !(isd:ied+1,jsd:jed  )
         integer, intent(IN):: npx,npy
         logical, OPTIONAL, intent(IN) :: DGRID, AGRID, CGRID, VECTOR
         integer :: i,j

         real :: mySign

         mySign = 1.0
         if (present(VECTOR)) then
            if (VECTOR) mySign = -1.0
         endif

         if (present(DGRID)) then
            call fill_corners_dgrid(x, y, npx, npy, mySign)
         elseif (present(CGRID)) then
            call fill_corners_cgrid(x, y, npx, npy, mySign)
         elseif (present(AGRID)) then
            call fill_corners_agrid(x, y, npx, npy, mySign)
         else
            call fill_corners_agrid(x, y, npx, npy, mySign)
         endif

      end subroutine fill_corners_xy_2d
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !!
!     fill_corners_xy_3d
      subroutine fill_corners_xy_3d(x, y, npx, npy, npz, DGRID, AGRID, CGRID, VECTOR)
         real, DIMENSION(isd:,jsd:,:), intent(INOUT):: x !(isd:ied  ,jsd:jed+1)
         real, DIMENSION(isd:,jsd:,:), intent(INOUT):: y !(isd:ied+1,jsd:jed  )
         integer, intent(IN):: npx,npy,npz
         logical, OPTIONAL, intent(IN) :: DGRID, AGRID, CGRID, VECTOR
         integer :: i,j,k

         real :: mySign

         mySign = 1.0
         if (present(VECTOR)) then
            if (VECTOR) mySign = -1.0
         endif

         if (present(DGRID)) then
            do k=1,npz
               call fill_corners_dgrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         elseif (present(CGRID)) then
            do k=1,npz
               call fill_corners_cgrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         elseif (present(AGRID)) then
            do k=1,npz
               call fill_corners_agrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         else
            do k=1,npz
               call fill_corners_agrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         endif

      end subroutine fill_corners_xy_3d
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
      subroutine fill_corners_dgrid(x, y, npx, npy, mySign)
         real, DIMENSION(isd:,jsd:), intent(INOUT):: x
         real, DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real, intent(IN) :: mySign 
         integer :: i,j

               do j=1,ng
                  do i=1,ng
                   !   if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) =        y(j+1  ,1-i    )  !SW Corner 
                   !   if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) = mySign*y(j+1  ,npy-1+i)  !NW Corner
                   !   if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) = mySign*y(npx-j,1-i    )  !SE Corner
                   !   if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) =        y(npx-j,npy-1+i)  !NE Corner
                      if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) = mySign*y(1-j  ,i    )  !SW Corner 
                      if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) =        y(1-j  ,npy-i)  !NW Corner
                      if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) =        y(npx+j,i    )  !SE Corner
                      if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) = mySign*y(npx+j,npy-i)  !NE Corner
                  enddo
               enddo
               do j=1,ng
                  do i=1,ng
                   !  if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) =        x(1-j    ,i+1  )  !SW Corner 
                   !  if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) = mySign*x(1-j    ,npy-i)  !NW Corner
                   !  if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) = mySign*x(npx-1+j,i+1  )  !SE Corner
                   !  if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) =        x(npx-1+j,npy-i)  !NE Corner
                     if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) = mySign*x(j      ,1-i  )  !SW Corner 
                     if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) =        x(j      ,npy+i)  !NW Corner
                     if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) =        x(npx-j  ,1-i  )  !SE Corner
                     if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) = mySign*x(npx-j  ,npy+i)  !NE Corner
                  enddo
               enddo

      end subroutine fill_corners_dgrid
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
      subroutine fill_corners_cgrid(x, y, npx, npy, mySign)
         real, DIMENSION(isd:,jsd:), intent(INOUT):: x
         real, DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real, intent(IN) :: mySign
         integer :: i,j

                  do j=1,ng
                     do i=1,ng
                        if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j    ) =        y(j      ,1-i  )  !SW Corner 
                        if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy-1+j) = mySign*y(j      ,npy+i)  !NW Corner
                        if ((ie+1==npx) .and. (js  ==  1)) x(npx+i  ,1-j    ) = mySign*y(npx-j  ,1-i  )  !SE Corner
                        if ((ie+1==npx) .and. (je+1==npy)) x(npx+i  ,npy-1+j) =        y(npx-j  ,npy+i)  !NE Corner
                     enddo
                  enddo
                  do j=1,ng
                     do i=1,ng
                        if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j  ) =        x(1-j  ,i    )  !SW Corner 
                        if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy+j) = mySign*x(1-j  ,npy-i)  !NW Corner
                        if ((ie+1==npx) .and. (js  ==  1)) y(npx-1+i,1-j  ) = mySign*x(npx+j,i    )  !SE Corner
                        if ((ie+1==npx) .and. (je+1==npy)) y(npx-1+i,npy+j) =        x(npx+j,npy-i)  !NE Corner
                     enddo
                  enddo
      
      end subroutine fill_corners_cgrid
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
      subroutine fill_corners_agrid(x, y, npx, npy, mySign)
         real, DIMENSION(isd:,jsd:), intent(INOUT):: x
         real, DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real, intent(IN) :: mySign
         integer :: i,j

                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) x(1-i    ,1-j    ) = mySign*y(1-j    ,i        )  !SW Corner
                       if ((is==    1) .and. (je==npy-1)) x(1-i    ,npy-1+j) =        y(1-j    ,npy-1-i+1)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) x(npx-1+i,1-j    ) =        y(npx-1+j,i        )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) x(npx-1+i,npy-1+j) = mySign*y(npx-1+j,npy-1-i+1)  !NE Corner
                    enddo
                 enddo
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) y(1-j    ,1-i    ) = mySign*x(i        ,1-j    )  !SW Corner
                       if ((is==    1) .and. (je==npy-1)) y(1-j    ,npy-1+i) =        x(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) y(npx-1+j,1-i    ) =        x(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) y(npx-1+j,npy-1+i) = mySign*x(npx-1-i+1,npy-1+j)  !NE Corner
                    enddo
                 enddo

      end subroutine fill_corners_agrid
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!     
!     mp_corner_comm :: Point-based MPI communcation routine for Cubed-Sphere
!                       ghosted corner point on B-Grid 
!                       this routine sends 24 16-byte messages 
!     
      subroutine mp_corner_comm(q, npx, npy)
         integer, intent(IN)  :: npx,npy
         real  , intent(INOUT):: q(isd:ied+1,jsd:jed+1)

         real :: qsend(24)
         real :: send_tag, recv_tag
         integer :: sqest(24), rqest(24)
         integer :: Stats(24*MPI_STATUS_SIZE)
         integer :: nsend, nrecv, nread
         integer :: dest_gid, src_gid
         integer :: n

         qsend = 1.e25
         nsend=0
         nrecv=0

         if ( mod(tile,2) == 0 ) then
! Even Face LL and UR pairs 6 2-way
            if ( (is==1) .and. (js==1) ) then
               nsend=nsend+1
               qsend(nsend) = q(is,js+1)
               send_tag = 300+tile
               dest_gid = (tile-2)*npes_x*npes_y - 1
               if (dest_gid < 0) dest_gid=npes+dest_gid
               recv_tag = 100+(tile-2)
               if (tile==2) recv_tag = 100+(ntiles)
               src_gid  = (tile-3)*npes_x*npes_y
               src_gid  = src_gid + npes_x*(npes_y-1) + npes_x - 1
               if (src_gid < 0) src_gid=npes+src_gid
               if (npes>6) then
                  call MPI_SENDRECV( qsend(nsend), 1, MPI_DOUBLE_PRECISION, &
                                     dest_gid, send_tag, &
                                     q(is-1,js), 1, MPI_DOUBLE_PRECISION, &
                                     src_gid, recv_tag, &
                                     commglobal, Stats, ierror )
                  nsend=nsend-1
               else 
                  call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &
                                  send_tag, commglobal, sqest(nsend), ierror )
                  nrecv=nrecv+1
                  call MPI_IRECV( q(is-1,js), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                                  recv_tag, commglobal, rqest(nrecv), ierror )
               endif
            endif
            if ( (ie==npx-1) .and. (je==npy-1) ) then
               nsend=nsend+1
               qsend(nsend) = q(ie,je+1)
               send_tag = 100+tile
               dest_gid = (tile+1)*npes_x*npes_y
               if (dest_gid+1 > npes) dest_gid=dest_gid-npes
               recv_tag = 300+(tile+2)
               if (tile==6) recv_tag = 300+2
               src_gid  = (tile+1)*npes_x*npes_y
               if (src_gid+1 > npes) src_gid=src_gid-npes
               if (npes>6) then
                  call MPI_SENDRECV( qsend(nsend), 1, MPI_DOUBLE_PRECISION, &
                                     dest_gid, send_tag, &
                                     q(ie+2,je+1), 1, MPI_DOUBLE_PRECISION, &
                                     src_gid, recv_tag, &
                                     commglobal, Stats, ierror )
                  nsend=nsend-1
               else
                  call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &
                                  send_tag, commglobal, sqest(nsend), ierror )
                  nrecv=nrecv+1
                  call MPI_IRECV( q(ie+2,je+1), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                                  recv_tag, commglobal, rqest(nrecv), ierror )
               endif
            endif
! wait for comm to complete
            if (npes==6) then
               if (nsend>0) then
                  call MPI_WAITALL(nsend, sqest, Stats, ierror)
!!$                  do n=1,nsend
!!$                     call MPI_REQUEST_FREE( sqest(n), ierror )
!!$                  enddo
               endif
               if (nrecv>0) then
                  call MPI_WAITALL(nrecv, rqest, Stats, ierror)
!!$                  do n=1,nrecv
!!$                     call MPI_REQUEST_FREE( rqest(n), ierror )
!!$                  enddo
               endif
               nsend=0 ; nrecv=0
            endif

! Even Face LR 1 pair ; 1 1-way
            if ( (tile==2) .and. (ie==npx-1) .and. (js==1) ) then
               nsend=nsend+1
               qsend(nsend) = q(ie,js)
               send_tag = 200+tile
               dest_gid = (tile+1)*npes_x*npes_y + npes_x-1
               recv_tag = 200+(tile+2)
               src_gid  = dest_gid
               if (npes>6) then
                  call MPI_SENDRECV( qsend(nsend), 1, MPI_DOUBLE_PRECISION, &
                                     dest_gid, send_tag, &
                                     q(ie+2,js), 1, MPI_DOUBLE_PRECISION, &
                                     src_gid, recv_tag, &
                                     commglobal, Stats, ierror )
                  nsend=nsend-1
               else
                  call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &
                                  send_tag, commglobal, sqest(nsend), ierror )
                  nrecv=nrecv+1
                  call MPI_IRECV( q(ie+2,js), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                                  recv_tag, commglobal, rqest(nrecv), ierror )
               endif
            endif
            if ( (tile==4) .and. (ie==npx-1) .and. (js==1) ) then 
               nsend=nsend+1
               qsend(nsend) = q(ie+1,js+1)
               send_tag = 200+tile
               dest_gid = (tile-3)*npes_x*npes_y + npes_x-1
               recv_tag = 200+(tile-2)
               src_gid  = dest_gid
               if (npes>6) then
                  call MPI_SENDRECV( qsend(nsend), 1, MPI_DOUBLE_PRECISION, &     
                                     dest_gid, send_tag, &
                                     q(ie+2,js), 1, MPI_DOUBLE_PRECISION, &
                                     src_gid, recv_tag, &          
                                     commglobal, Stats, ierror )
                  nsend=nsend-1
               else
                  call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &        
                                  send_tag, commglobal, sqest(nsend), ierror )
                  nrecv=nrecv+1
                  call MPI_IRECV( q(ie+2,js), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                                  recv_tag, commglobal, rqest(nrecv), ierror )
               endif
               nsend=nsend+1
               qsend(nsend) = q(ie,js)
               send_tag = 200+tile
               dest_gid = (tile+1)*npes_x*npes_y + npes_x-1
               call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &
                               send_tag, commglobal, sqest(nsend), ierror )
            endif
            if ( (tile==6) .and. (ie==npx-1) .and. (js==1) ) then
               recv_tag = 200+(tile-2)
               src_gid  = (tile-3)*npes_x*npes_y + npes_x-1
               nrecv=nrecv+1
               call MPI_IRECV( q(ie+2,js), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                               recv_tag, commglobal, rqest(nrecv), ierror )
            endif

! wait for comm to complete 
            if (npes==6) then
               if (nsend>0) then
                  call MPI_WAITALL(nsend, sqest, Stats, ierror)
!!$                  do n=1,nsend
!!$                     call MPI_REQUEST_FREE( sqest(n), ierror )
!!$                  enddo
               endif
               if (nrecv>0) then
                  call MPI_WAITALL(nrecv, rqest, Stats, ierror)
!!$                  do n=1,nrecv
!!$                     call MPI_REQUEST_FREE( rqest(n), ierror )
!!$                  enddo
               endif
               nsend=0 ; nrecv=0
            endif

! Send to Odd face LR 3 1-way
            if ( (is==1) .and. (js==1) ) then
               nsend=nsend+1
               qsend(nsend) = q(is+1,js)
               send_tag = 200+tile
               dest_gid = (tile-2)*npes_x*npes_y + npes_x-1
               call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &
                               send_tag, commglobal, sqest(nsend), ierror )
            endif

! Receive Even Face UL 3 1-way
            if ( (is==1) .and. (je==npy-1) ) then
               recv_tag = 400+(tile-1)
               src_gid  = (tile-2)*npes_x*npes_y + npes_x*(npes_y-1) + npes_x-1
               nrecv=nrecv+1
               call MPI_IRECV( q(is-1,je+1), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                               recv_tag, commglobal, rqest(nrecv), ierror )
            endif

         else

! Odd Face LL and UR pairs 6 2-way
            if ( (is==1) .and. (js==1) ) then
               nsend=nsend+1
               qsend(nsend) = q(is+1,js)
               send_tag = 300+tile
               dest_gid = (tile-2)*npes_x*npes_y - 1
               if (dest_gid < 0) dest_gid=npes+dest_gid
               recv_tag = 100+(tile-2)
               if (tile==1) recv_tag = 100+(ntiles-tile)
               src_gid  = (tile-3)*npes_x*npes_y
               src_gid  = src_gid + npes_x*(npes_y-1) + npes_x - 1
               if (src_gid < 0) src_gid=npes+src_gid
               if (npes>6) then
                  call MPI_SENDRECV( qsend(nsend), 1, MPI_DOUBLE_PRECISION, &
                                     dest_gid, send_tag, &
                                     q(is-1,js), 1, MPI_DOUBLE_PRECISION, &
                                     src_gid, recv_tag, &             
                                     commglobal, Stats, ierror )
                  nsend=nsend-1
               else 
                  call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &
                                  send_tag, commglobal, sqest(nsend), ierror )
                  nrecv=nrecv+1
                  call MPI_IRECV( q(is-1,js), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                                  recv_tag, commglobal, rqest(nrecv), ierror )
               endif
            endif
            if ( (ie==npx-1) .and. (je==npy-1) ) then
               nsend=nsend+1
               qsend(nsend) = q(ie+1,je)
               send_tag = 100+tile
               dest_gid = (tile+1)*npes_x*npes_y
               if (dest_gid+1 > npes) dest_gid=dest_gid-npes
               recv_tag = 300+(tile+2)
               if (tile==5) recv_tag = 300+1
               src_gid  = (tile+1)*npes_x*npes_y
               if (src_gid+1 > npes) src_gid=src_gid-npes
               if (npes>6) then
                  call MPI_SENDRECV( qsend(nsend), 1, MPI_DOUBLE_PRECISION, &      
                                     dest_gid, send_tag, &
                                     q(ie+2,je+1), 1, MPI_DOUBLE_PRECISION, &
                                     src_gid, recv_tag, &             
                                     commglobal, Stats, ierror )
                  nsend=nsend-1
               else 
                  call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &
                                  send_tag, commglobal, sqest(nsend), ierror )
                  nrecv=nrecv+1
                  call MPI_IRECV( q(ie+2,je+1), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                                  recv_tag, commglobal, rqest(nrecv), ierror )
               endif
            endif
! wait for comm to complete 
            if (npes==6) then
               if (nsend>0) then
                  call MPI_WAITALL(nsend, sqest, Stats, ierror)
!!$                  do n=1,nsend
!!$                     call MPI_REQUEST_FREE( sqest(n), ierror )
!!$                  enddo
               endif
               if (nrecv>0) then
                  call MPI_WAITALL(nrecv, rqest, Stats, ierror)
!!$                  do n=1,nrecv
!!$                     call MPI_REQUEST_FREE( rqest(n), ierror )
!!$                  enddo
               endif
               nsend=0 ; nrecv=0
            endif
            
! Odd Face UL 1 pair ; 1 1-way
            if ( (tile==1) .and. (is==1) .and. (je==npy-1) ) then
               nsend=nsend+1
               qsend(nsend) = q(is,je)
               send_tag = 400+tile
               dest_gid = (tile+1)*npes_x*npes_y + npes_x*(npes_y-1)
               recv_tag = 400+(tile+2)
               src_gid  = dest_gid
               if (npes>6) then
                  call MPI_SENDRECV( qsend(nsend), 1, MPI_DOUBLE_PRECISION, &
                                     dest_gid, send_tag, &
                                     q(is-1,je+1), 1, MPI_DOUBLE_PRECISION, &
                                     src_gid, recv_tag, &
                                     commglobal, Stats, ierror )
                  nsend=nsend-1
               else
                  call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &
                                  send_tag, commglobal, sqest(nsend), ierror )
                  nrecv=nrecv+1
                  call MPI_IRECV( q(is-1,je+1), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                                  recv_tag, commglobal, rqest(nrecv), ierror )
               endif
            endif
            if ( (tile==3) .and. (is==1) .and. (je==npy-1) ) then
               nsend=nsend+1
               qsend(nsend) = q(is+1,je+1)
               send_tag = 400+tile
               dest_gid = npes_x*(npes_y-1)
               recv_tag = 400+(tile-2)           
               src_gid  = dest_gid
               if (npes>6) then
                  call MPI_SENDRECV( qsend(nsend), 1, MPI_DOUBLE_PRECISION, &     
                                     dest_gid, send_tag, &
                                     q(is-1,je+1), 1, MPI_DOUBLE_PRECISION, &
                                     src_gid, recv_tag, &          
                                     commglobal, Stats, ierror )
                  nsend=nsend-1             
               else
                  call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &        
                                  send_tag, commglobal, sqest(nsend), ierror )
                  nrecv=nrecv+1
                  call MPI_IRECV( q(is-1,je+1), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                                  recv_tag, commglobal, rqest(nrecv), ierror )
               endif            
               nsend=nsend+1
               qsend(nsend) = q(is,je)
               send_tag = 400+tile
               dest_gid = (tile+1)*npes_x*npes_y + npes_x*(npes_y-1)
               call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &
                               send_tag, commglobal, sqest(nsend), ierror )
            endif
            if ( (tile==5) .and. (is==1) .and. (je==npy-1) ) then
               recv_tag = 400+(tile-2)
               src_gid  = (tile-3)*npes_x*npes_y + npes_x*(npes_y-1)
               nrecv=nrecv+1 
               call MPI_IRECV( q(is-1,je+1), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                               recv_tag, commglobal, rqest(nrecv), ierror ) 
            endif

! wait for comm to complete 
            if (npes==6) then
               if (nsend>0) then
                  call MPI_WAITALL(nsend, sqest, Stats, ierror)
!!$                  do n=1,nsend
!!$                     call MPI_REQUEST_FREE( sqest(n), ierror )
!!$                  enddo
               endif
               if (nrecv>0) then
                  call MPI_WAITALL(nrecv, rqest, Stats, ierror)
!!$                  do n=1,nrecv
!!$                     call MPI_REQUEST_FREE( rqest(n), ierror )
!!$                  enddo
               endif
               nsend=0 ; nrecv=0
            endif

! Send to Even face UL 3 1-way 
            if ( (ie==npx-1) .and. (je==npy-1) ) then
               nsend=nsend+1
               qsend(nsend) = q(ie,je+1)
               send_tag = 400+tile
               dest_gid = tile*npes_x*npes_y + npes_x*(npes_y-1)
               call MPI_ISEND( qsend(nsend), 1, MPI_DOUBLE_PRECISION, dest_gid, &
                               send_tag, commglobal, sqest(nsend), ierror )
            endif

! Receive Odd Face LR 3 1-way
            if ( (ie==npx-1) .and. (js==1) ) then
               recv_tag = 200+(tile+1)
               src_gid  = (tile-1)*npes_x*npes_y + npes_x*npes_y 
               nrecv=nrecv+1
               call MPI_IRECV( q(ie+2,js), 1, MPI_DOUBLE_PRECISION, src_gid,  &
                               recv_tag, commglobal, rqest(nrecv), ierror )
            endif

         endif

! wait for comm to complete
         if (nsend>0) then
            call MPI_WAITALL(nsend, sqest, Stats, ierror)
!!$            do n=1,nsend
!!$               call MPI_REQUEST_FREE( sqest(n), ierror )
!!$            enddo
         endif
         if (nrecv>0) then
            call MPI_WAITALL(nrecv, rqest, Stats, ierror)
!!$            do n=1,nrecv
!!$               call MPI_REQUEST_FREE( rqest(n), ierror )
!!$            enddo
         endif

      end subroutine mp_corner_comm
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!       
!     mp_gather_4d_r4 :: Call SPMD Gather 
!     
      subroutine mp_gather_4d_r4(q, i1,i2, j1,j2, idim, jdim, kdim, ldim)
         integer, intent(IN)  :: i1,i2, j1,j2
         integer, intent(IN)  :: idim, jdim, kdim, ldim
         real(kind=4), intent(INOUT):: q(idim,jdim,kdim,ldim)
         integer :: i,j,k,l,n,icnt 
         integer :: Lsize, Lsize_buf(1)
         integer :: Gsize
         integer :: LsizeS(npes), Ldispl(npes), cnts(npes)
         integer :: Ldims(5), Gdims(5*npes)
         real(kind=4), allocatable, dimension(:) :: larr, garr
        
         Ldims(1) = i1
         Ldims(2) = i2
         Ldims(3) = j1
         Ldims(4) = j2
         Ldims(5) = tile 
         do l=1,npes
            cnts(l) = 5
            Ldispl(l) = 5*(l-1)
         enddo 
         call mpp_gather(Ldims, Gdims)
!         call MPI_GATHERV(Ldims, 5, MPI_INTEGER, Gdims, cnts, Ldispl, MPI_INTEGER, masterproc, commglobal, ierror)
      
         Lsize = ( (i2 - i1 + 1) * (j2 - j1 + 1) ) * kdim
         do l=1,npes
            cnts(l) = 1
            Ldispl(l) = l-1
         enddo 
         LsizeS(:)=1
         Lsize_buf(1) = Lsize
         call mpp_gather(Lsize_buf, LsizeS)
!         call MPI_GATHERV(Lsize, 1, MPI_INTEGER, LsizeS, cnts, Ldispl, MPI_INTEGER, masterproc, commglobal, ierror)

         allocate ( larr(Lsize) )
         icnt = 1
         do k=1,kdim
            do j=j1,j2
               do i=i1,i2
                  larr(icnt) = q(i,j,k,tile)
                  icnt=icnt+1
               enddo
            enddo
         enddo
         Ldispl(1) = 0.0
!         call mp_bcst(LsizeS(1))
         call mpp_broadcast(LsizeS, npes, masterproc)
         Gsize = LsizeS(1)
         do l=2,npes
!            call mp_bcst(LsizeS(l))
            Ldispl(l) = Ldispl(l-1) + LsizeS(l-1)
            Gsize = Gsize + LsizeS(l)
         enddo
         allocate ( garr(Gsize) )

         call mpp_gather(larr, Lsize, garr, LsizeS)
!         call MPI_GATHERV(larr, Lsize, MPI_REAL, garr, LsizeS, Ldispl, MPI_REAL, masterproc, commglobal, ierror)

         if (gid==masterproc) then
            do n=2,npes
               icnt=1
               do l=Gdims( (n-1)*5 + 5 ), Gdims( (n-1)*5 + 5 )
                  do k=1,kdim
                     do j=Gdims( (n-1)*5 + 3 ), Gdims( (n-1)*5 + 4 )
                        do i=Gdims( (n-1)*5 + 1 ), Gdims( (n-1)*5 + 2 )
                           q(i,j,k,l) = garr(Ldispl(n)+icnt)
                           icnt=icnt+1
                        enddo
                    enddo
                 enddo
               enddo
            enddo
         endif
         deallocate( larr )
         deallocate( garr )

      end subroutine mp_gather_4d_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_gather_3d_r4 :: Call SPMD Gather 
!
      subroutine mp_gather_3d_r4(q, i1,i2, j1,j2, idim, jdim, ldim)
         integer, intent(IN)  :: i1,i2, j1,j2
         integer, intent(IN)  :: idim, jdim, ldim
         real(kind=4), intent(INOUT):: q(idim,jdim,ldim)
         integer :: i,j,l,n,icnt 
         integer :: Lsize, Lsize_buf(1)
         integer :: Gsize
         integer :: LsizeS(npes), Ldispl(npes), cnts(npes)
         integer :: Ldims(5), Gdims(5*npes)
         real(kind=4), allocatable, dimension(:) :: larr, garr 

         Ldims(1) = i1
         Ldims(2) = i2
         Ldims(3) = j1
         Ldims(4) = j2
         Ldims(5) = tile
         do l=1,npes
            cnts(l) = 5
            Ldispl(l) = 5*(l-1)
         enddo
!         call MPI_GATHERV(Ldims, 5, MPI_INTEGER, Gdims, cnts, Ldispl, MPI_INTEGER, masterproc, commglobal, ierror)
         call mpp_gather(Ldims, Gdims)

         Lsize = ( (i2 - i1 + 1) * (j2 - j1 + 1) )
         do l=1,npes
            cnts(l) = 1
            Ldispl(l) = l-1
         enddo 
         LsizeS(:)=1
         Lsize_buf(1) = Lsize
         call mpp_gather(Lsize_buf, LsizeS)
!         call MPI_GATHERV(Lsize, 1, MPI_INTEGER, LsizeS, cnts, Ldispl, MPI_INTEGER, masterproc, commglobal, ierror)

         allocate ( larr(Lsize) )
         icnt = 1
         do j=j1,j2
            do i=i1,i2
               larr(icnt) = q(i,j,tile)  
               icnt=icnt+1
            enddo
         enddo
         Ldispl(1) = 0.0
!         call mp_bcst(LsizeS(1))
         call mpp_broadcast(LsizeS, npes, masterproc)
         Gsize = LsizeS(1)
         do l=2,npes
!            call mp_bcst(LsizeS(l))
            Ldispl(l) = Ldispl(l-1) + LsizeS(l-1)
            Gsize = Gsize + LsizeS(l)
         enddo
         allocate ( garr(Gsize) )
         call mpp_gather(larr, Lsize, garr, LsizeS)
!         call MPI_GATHERV(larr, Lsize, MPI_REAL, garr, LsizeS, Ldispl, MPI_REAL, masterproc, commglobal, ierror)
         if (gid==masterproc) then
            do n=2,npes
               icnt=1
               do l=Gdims( (n-1)*5 + 5 ), Gdims( (n-1)*5 + 5 )
                  do j=Gdims( (n-1)*5 + 3 ), Gdims( (n-1)*5 + 4 ) 
                     do i=Gdims( (n-1)*5 + 1 ), Gdims( (n-1)*5 + 2 )
                        q(i,j,l) = garr(Ldispl(n)+icnt)
                        icnt=icnt+1
                     enddo
                 enddo
               enddo
            enddo
         endif
         deallocate( larr )
         deallocate( garr )

      end subroutine mp_gather_3d_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_gather_3d_r8 :: Call SPMD Gather 
!
      subroutine mp_gather_3d_r8(q, i1,i2, j1,j2, idim, jdim, ldim)
         integer, intent(IN)  :: i1,i2, j1,j2
         integer, intent(IN)  :: idim, jdim, ldim
         real(kind=8),  intent(INOUT):: q(idim,jdim,ldim)
         integer :: i,j,l,n,icnt
         integer :: Lsize, Lsize_buf(1)
         integer :: Gsize
         integer :: LsizeS(npes), Ldispl(npes), cnts(npes)
         integer :: Ldims(5), Gdims(5*npes)
         real,   allocatable, dimension(:) :: larr, garr

         Ldims(1) = i1
         Ldims(2) = i2
         Ldims(3) = j1
         Ldims(4) = j2
         Ldims(5) = tile
         do l=1,npes
            cnts(l) = 5
            Ldispl(l) = 5*(l-1)
         enddo
!         call MPI_GATHER(Ldims, 5, MPI_INTEGER, Gdims, cnts, MPI_INTEGER, masterproc, commglobal, ierror)
         call mpp_gather(Ldims, Gdims)
         Lsize = ( (i2 - i1 + 1) * (j2 - j1 + 1) )
         do l=1,npes
            cnts(l) = 1
            Ldispl(l) = l-1
         enddo
         LsizeS(:)=0.

!         call MPI_GATHERV(Lsize, 1, MPI_INTEGER, LsizeS, cnts, Ldispl, MPI_INTEGER, masterproc, commglobal, ierror)
         Lsize_buf(1) = Lsize
         call mpp_gather(Lsize_buf, LsizeS)

         allocate ( larr(Lsize) )
         icnt = 1
         do j=j1,j2
            do i=i1,i2
               larr(icnt) = q(i,j,tile)
               icnt=icnt+1
            enddo
         enddo
         Ldispl(1) = 0.0
         call mpp_broadcast(LsizeS, npes, masterproc)
!         call mp_bcst(LsizeS(1))
         Gsize = LsizeS(1)
         do l=2,npes
!            call mp_bcst(LsizeS(l))
            Ldispl(l) = Ldispl(l-1) + LsizeS(l-1)
            Gsize = Gsize + LsizeS(l)
         enddo

         allocate ( garr(Gsize) )
         call mpp_gather(larr, Lsize, garr, LsizeS)
!         call MPI_GATHERV(larr, Lsize, MPI_DOUBLE_PRECISION, garr, LsizeS, Ldispl, MPI_DOUBLE_PRECISION, masterproc, commglobal, ierror)
         if (gid==masterproc) then
            do n=2,npes
               icnt=1
               do l=Gdims( (n-1)*5 + 5 ), Gdims( (n-1)*5 + 5 )
                  do j=Gdims( (n-1)*5 + 3 ), Gdims( (n-1)*5 + 4 )
                     do i=Gdims( (n-1)*5 + 1 ), Gdims( (n-1)*5 + 2 )
                        q(i,j,l) = garr(Ldispl(n)+icnt)
                        icnt=icnt+1
                     enddo
                 enddo
               enddo
            enddo
         endif
         deallocate( larr )
         deallocate( garr )

      end subroutine mp_gather_3d_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_i4 :: Call SPMD broadcast 
!
      subroutine mp_bcst_i4(q)
         integer, intent(INOUT)  :: q

         call MPI_BCAST(q, 1, MPI_INTEGER, masterproc, commglobal, ierror)

      end subroutine mp_bcst_i4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_r8 :: Call SPMD broadcast 
!
      subroutine mp_bcst_r8(q)
         real, intent(INOUT)  :: q

         call MPI_BCAST(q, 1, MPI_DOUBLE_PRECISION, masterproc, commglobal, ierror)

      end subroutine mp_bcst_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_3d_r8 :: Call SPMD broadcast 
!
      subroutine mp_bcst_3d_r8(q, idim, jdim, kdim)
         integer, intent(IN)  :: idim, jdim, kdim
         real, intent(INOUT)  :: q(idim,jdim,kdim)

         call MPI_BCAST(q, idim*jdim*kdim, MPI_DOUBLE_PRECISION, masterproc, commglobal, ierror)

      end subroutine mp_bcst_3d_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!       
!     mp_bcst_4d_r8 :: Call SPMD broadcast 
!
      subroutine mp_bcst_4d_r8(q, idim, jdim, kdim, ldim)
         integer, intent(IN)  :: idim, jdim, kdim, ldim
         real, intent(INOUT)  :: q(idim,jdim,kdim,ldim)

         call MPI_BCAST(q, idim*jdim*kdim*ldim, MPI_DOUBLE_PRECISION, masterproc, commglobal, ierror)
        
      end subroutine mp_bcst_4d_r8
!     
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!       
!     mp_reduce_max_r8_1d :: Call SPMD REDUCE_MAX 
!
      subroutine mp_reduce_max_r8_1d(mymax,npts)
         integer, intent(IN)  :: npts
         real, intent(INOUT)  :: mymax(npts)
        
         real :: gmax(npts)
        
         call MPI_ALLREDUCE( mymax, gmax, npts, MPI_DOUBLE_PRECISION, MPI_MAX, &
                             commglobal, ierror )
      
         mymax = gmax
        
      end subroutine mp_reduce_max_r8_1d
!     
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_reduce_max_r8 :: Call SPMD REDUCE_MAX 
!
      subroutine mp_reduce_max_r8(mymax)
         real, intent(INOUT)  :: mymax

         real :: gmax

         call MPI_ALLREDUCE( mymax, gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                             commglobal, ierror )

         mymax = gmax

      end subroutine mp_reduce_max_r8

      subroutine mp_reduce_min(mymin)
         real, intent(INOUT)  :: mymin

         real :: gmin

         call MPI_ALLREDUCE( mymin, gmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                             commglobal, ierror )

         mymin = gmin

      end subroutine mp_reduce_min
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_4d_i4 :: Call SPMD REDUCE_MAX 
!
      subroutine mp_reduce_max_i4(mymax)
         integer, intent(INOUT)  :: mymax

         integer :: gmax

         call MPI_ALLREDUCE( mymax, gmax, 1, MPI_INTEGER, MPI_MAX, &
                             commglobal, ierror )

         mymax = gmax

      end subroutine mp_reduce_max_i4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_reduce_sum_r8 :: Call SPMD REDUCE_SUM 
!
      subroutine mp_reduce_sum_r8(mysum)
         real, intent(INOUT)  :: mysum

         real :: gsum

         call MPI_ALLREDUCE( mysum, gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                             commglobal, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_reduce_sum_r8_1d :: Call SPMD REDUCE_SUM 
!
      subroutine mp_reduce_sum_r8_1d(mysum, sum1d, npts)
         integer, intent(in)  :: npts
         real, intent(in)     :: sum1d(npts)
         real, intent(INOUT)  :: mysum

         real :: gsum
         integer :: i

         mysum = 0.0
         do i=1,npts
            mysum = mysum + sum1d(i)
         enddo 

         call MPI_ALLREDUCE( mysum, gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                             commglobal, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_r8_1d
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

#else
      implicit none
      private
      integer :: masterproc = 0
      integer :: gid = 0
      integer, parameter:: ng    = 3     ! Number of ghost zones required
      public gid, masterproc, ng
#endif

      end module fv_mp_mod
!-------------------------------------------------------------------------------
