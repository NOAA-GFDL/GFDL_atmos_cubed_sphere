module fv_tracer2d_mod
   use tp_core_mod,       only: fv_tp_2d, copy_corners
   use fv_mp_mod,         only: mp_reduce_max
   use fv_mp_mod,         only: ng, mp_gather, is_master
   use mpp_domains_mod,   only: mpp_start_update_domains, mpp_complete_update_domains, mpp_update_domains, CGRID_NE, domain2d
   use fv_timing_mod,     only: timing_on, timing_off
   use boundary_mod,      only: nested_grid_BC_apply, nested_grid_BC_apply_intT
   use fv_arrays_mod,     only: fv_grid_type, fv_nest_type, fv_atmos_type, fv_grid_bounds_type
   use mpp_mod,           only: mpp_error, FATAL, mpp_broadcast, mpp_send, mpp_recv, mpp_sum, mpp_max

implicit none
private

public :: tracer_2d, tracer_2d_nested, tracer_2d_1L

real, allocatable, dimension(:,:,:) :: nest_fx_west_accum, nest_fx_east_accum, nest_fx_south_accum, nest_fx_north_accum

!---- version number -----
   character(len=128) :: version = '$Id: fv_tracer2d.F90,v 20.0 2013/12/13 23:04:36 fms Exp $'
   character(len=128) :: tagname = '$Name: tikal_201409 $'

contains

!-----------------------------------------------------------------------
! !ROUTINE: Perform 2D horizontal-to-lagrangian transport
!-----------------------------------------------------------------------

subroutine tracer_2d_1L(q, dp0, mfx, mfy, cx, cy, gridstruct, neststruct, bd, domain, npx, npy, npz, nq, hord,  &
     q_split, k, q3, dt, id_divg, k_split)
      type(fv_grid_bounds_type), intent(IN) :: bd
      integer, intent(IN) :: npx, npy, npz
      integer, intent(IN) :: k
      integer, intent(IN) :: nq    ! number of tracers to be advected
      integer, intent(IN) :: hord
      integer, intent(IN) :: q_split, k_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt
      real   , intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,nq)       ! 2D Tracers
      real   , intent(INOUT) ::q3(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)   ! Tracers
      real   , intent(INOUT) :: dp0(bd%is:bd%ie,bd%js:bd%je)        ! DELP before dyn_core
      real   , intent(IN) :: mfx(bd%is:bd%ie+1,bd%js:bd%je)    ! Mass Flux X-Dir
      real   , intent(IN) :: mfy(bd%is:bd%ie  ,bd%js:bd%je+1)    ! Mass Flux Y-Dir
      real   , intent(IN) ::  cx(bd%is:bd%ie+1,bd%jsd:bd%jed)  ! Courant Number X-Dir
      real   , intent(IN) ::  cy(bd%isd:bd%ied,bd%js :bd%je +1)  ! Courant Number Y-Dir
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(fv_nest_type), intent(INOUT) :: neststruct
      type(domain2d), intent(INOUT) :: domain

! Local Arrays
      real :: mfx2(bd%is:bd%ie+1,bd%js:bd%je)
      real :: mfy2(bd%is:bd%ie  ,bd%js:bd%je+1)
      real ::  cx2(bd%is:bd%ie+1,bd%jsd:bd%jed)
      real ::  cy2(bd%isd:bd%ied,bd%js :bd%je +1)

      real :: dp1(bd%is:bd%ie,bd%js:bd%je)
      real :: dp2(bd%is:bd%ie,bd%js:bd%je)
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je )
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1)
      real :: ra_x(bd%is:bd%ie,bd%jsd:bd%jed)
      real :: ra_y(bd%isd:bd%ied,bd%js:bd%je)
      real :: xfx(bd%is:bd%ie+1,bd%jsd:bd%jed)
      real :: yfx(bd%isd:bd%ied,bd%js: bd%je+1)
      real :: cmax
      real :: frac, rdt
      integer :: nsplt
      integer :: i,j,it,iq
      integer :: i_pack

      real, pointer, dimension(:,:) :: area, rarea, sina_u, sina_v
      real, pointer, dimension(:,:) :: dxa, dya, dx, dy
      real, pointer, dimension(:,:,:) :: sin_sg

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

       area => gridstruct%area
      rarea => gridstruct%rarea

      sina_u => gridstruct%sina_u
      sina_v => gridstruct%sina_v

      sin_sg => gridstruct%sin_sg
      dxa    => gridstruct%dxa 
      dya    => gridstruct%dya 
      dx     => gridstruct%dx  
      dy     => gridstruct%dy  

      i_pack = mpp_start_update_domains(q, domain)

      do j=jsd,jed
         do i=is,ie+1
            if (cx(i,j) > 0.) then
                xfx(i,j) = cx(i,j)*dxa(i-1,j)*dy(i,j)*sin_sg(i-1,j,3)
            else
                xfx(i,j) = cx(i,j)*dxa(i,j)*dy(i,j)*sin_sg(i,j,1)
            endif
         enddo
      enddo

      do j=js,je+1
         do i=isd,ied
            if (cy(i,j) > 0.) then
                yfx(i,j) = cy(i,j)*dya(i,j-1)*dx(i,j)*sin_sg(i,j-1,4)
            else
                yfx(i,j) = cy(i,j)*dya(i,j)*dx(i,j)*sin_sg(i,j,2)
            endif
         enddo
      enddo


      if ( q_split==0 ) then
! Determine nsplt for tracer advection
         cmax = 0.
         do j=js,je
            do i=is,ie
               cmax = max(abs(cx(i,j))+(1.-sina_u(i,j)),     &
                          abs(cy(i,j))+(1.-sina_v(i,j)), cmax)
            enddo
         enddo
         call mp_reduce_max(cmax)
         nsplt = int(1.01 + cmax)
         if ( is_master() .and. nsplt > 5 )  write(*,*) k, 'Tracer_2d_split=', nsplt, cmax
      else
         nsplt = q_split
      endif

      frac  = 1. / real(nsplt)

          do j=jsd,jed
             do i=is,ie+1
                cx2(i,j) =  cx(i,j) * frac
                xfx(i,j) = xfx(i,j) * frac
             enddo
          enddo

          do j=js,je
             do i=is,ie+1
                mfx2(i,j) = mfx(i,j) * frac
             enddo
          enddo

          do j=js,je+1
             do i=isd,ied
                cy2(i,j) =  cy(i,j) * frac
               yfx(i,j) = yfx(i,j) * frac
             enddo
          enddo

          do j=js,je+1
             do i=is,ie
                mfy2(i,j) = mfy(i,j) * frac
             enddo
          enddo

      do j=jsd,jed
         do i=is,ie
            ra_x(i,j) = area(i,j) + xfx(i,j) - xfx(i+1,j)
         enddo
      enddo

      do j=js,je
         do i=isd,ied
            ra_y(i,j) = area(i,j) + yfx(i,j) - yfx(i,j+1)
         enddo
      enddo

      do j=js,je
         do i=is,ie
            dp1(i,j) = dp0(i,j)
         enddo
      enddo

! Start time split ...
      do it=1,nsplt

         do j=js,je
            do i=is,ie
               dp2(i,j) = dp1(i,j) + (mfx2(i,j) - mfx2(i+1,j) +  &
                          mfy2(i,j) - mfy2(i,j+1)) * rarea(i,j)
            enddo
         enddo

         call mpp_complete_update_domains(i_pack, q, domain)

         do iq=1,nq
            call fv_tp_2d( q(isd,jsd,iq), cx2, cy2, npx, npy, hord, fx, fy, &
                           xfx, yfx, gridstruct, bd, ra_x, ra_y, mfx=mfx2, mfy=mfy2 )
            if( it==nsplt ) then
            do j=js,je
               do i=is,ie
                  q3(i,j,k,iq) = (q(i,j,iq)*dp1(i,j) + (fx(i,j)-fx(i+1,j) + &
                                  fy(i,j)-fy(i,j+1))*rarea(i,j)) / dp2(i,j)
               enddo
            enddo
            else
            do j=js,je
               do i=is,ie
                  q(i,j,iq) = (q(i,j,iq)*dp1(i,j) + (fx(i,j)-fx(i+1,j) + &
                              fy(i,j)-fy(i,j+1))*rarea(i,j)) / dp2(i,j)
               enddo
            enddo
           endif

           !Apply nested-grid BCs
           if ( gridstruct%nested ) then

                 call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,iq), &
                      0, 0, npx, npy, bd, &
                      real(neststruct%tracer_nest_timestep)+real(nsplt*k_split), real(nsplt*k_split), &
                      var_east_t0=neststruct%q_BC(iq)%east_t0(:,:,k), &
                      var_west_t0=neststruct%q_BC(iq)%west_t0(:,:,k), &
                      var_north_t0=neststruct%q_BC(iq)%north_t0(:,:,k), &
                      var_south_t0=neststruct%q_BC(iq)%south_t0(:,:,k), &
                      var_east_t1=neststruct%q_BC(iq)%east_t1(:,:,k), &
                      var_west_t1=neststruct%q_BC(iq)%west_t1(:,:,k), &
                      var_north_t1=neststruct%q_BC(iq)%north_t1(:,:,k), &
                      var_south_t1=neststruct%q_BC(iq)%south_t1(:,:,k), &
                      bctype=neststruct%nestbctype, &
                      nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )
           end if

        enddo!q-loop

         if ( it/=nsplt ) then
              i_pack = mpp_start_update_domains(q, domain)
              do j=js,je
                 do i=is,ie
                    dp1(i,j) = dp2(i,j)
                 enddo
              enddo
         endif
     enddo  ! nsplt

     if ( gridstruct%nested .and. k == npz) then
        neststruct%tracer_nest_timestep = neststruct%tracer_nest_timestep + 1
     end if

     if ( id_divg > 0 ) then
         rdt = 1./(frac*dt)

         do j=js,je
            do i=is,ie
               dp0(i,j) = (xfx(i+1,j)-xfx(i,j) + yfx(i,j+1)-yfx(i,j))*rarea(i,j)*rdt
            enddo
         enddo
     endif

end subroutine tracer_2d_1L


subroutine tracer_2d(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz,   &
                     nq,  hord, q_split, dt, id_divg, q_pack, z_tracer, k_split)

      type(fv_grid_bounds_type), intent(IN) :: bd
      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    ! number of tracers to be advected
      integer, intent(IN) :: hord
      integer, intent(IN) :: q_split, k_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt
      logical, intent(IN) :: z_tracer
      integer, intent(inout) :: q_pack
      real   , intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)   ! Tracers
      real   , intent(INOUT) :: dp1(bd%is:bd%ie,bd%js:bd%je,npz)        ! DELP before dyn_core
      real   , intent(INOUT) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,  npz)    ! Mass Flux X-Dir
      real   , intent(INOUT) :: mfy(bd%is:bd%ie  ,bd%js:bd%je+1,npz)    ! Mass Flux Y-Dir
      real   , intent(INOUT) ::  cx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)  ! Courant Number X-Dir
      real   , intent(INOUT) ::  cy(bd%isd:bd%ied,bd%js :bd%je +1,npz)  ! Courant Number Y-Dir
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain

! Local Arrays
      real :: dp2(bd%is:bd%ie,bd%js:bd%je,npz)
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je )
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1)
      real :: ra_x(bd%is:bd%ie,bd%jsd:bd%jed)
      real :: ra_y(bd%isd:bd%ied,bd%js:bd%je)
      real :: xfx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)
      real :: yfx(bd%isd:bd%ied,bd%js: bd%je+1, npz)
      real :: cmax(npz)
      real :: cmax_t
      real :: c_global
      real :: frac, rdt
      integer :: nsplt
      integer :: i,j,k,it,iq

      real, pointer, dimension(:,:) :: area, rarea
      real, pointer, dimension(:,:,:) :: sin_sg
      real, pointer, dimension(:,:) :: dxa, dya, dx, dy

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

       area => gridstruct%area
      rarea => gridstruct%rarea

      sin_sg => gridstruct%sin_sg
      dxa    => gridstruct%dxa 
      dya    => gridstruct%dya 
      dx     => gridstruct%dx  
      dy     => gridstruct%dy  

!$omp parallel do default(shared)
      do k=1,npz
         do j=jsd,jed
            do i=is,ie+1
               if (cx(i,j,k) > 0.) then
                  xfx(i,j,k) = cx(i,j,k)*dxa(i-1,j)*dy(i,j)*sin_sg(i-1,j,3)
               else
                  xfx(i,j,k) = cx(i,j,k)*dxa(i,j)*dy(i,j)*sin_sg(i,j,1)
               endif
            enddo
         enddo
         do j=js,je+1
            do i=isd,ied
               if (cy(i,j,k) > 0.) then
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j-1)*dx(i,j)*sin_sg(i,j-1,4)
               else
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j)*dx(i,j)*sin_sg(i,j,2)
               endif
            enddo
         enddo
      enddo

!--------------------------------------------------------------------------------
  if ( q_split == 0 ) then
! Determine nsplt
!$omp parallel do default(shared) private(cmax_t )
      do k=1,npz
         cmax(k) = 0.
         if ( k < 4 ) then
! Top layers: C < max( abs(c_x), abs(c_y) )
            do j=js,je
               do i=is,ie
                  cmax_t  = max( abs(cx(i,j,k)), abs(cy(i,j,k)) )
                  cmax(k) = max( cmax_t, cmax(k) )
               enddo
            enddo
         else
            do j=js,je
               do i=is,ie
                  cmax_t  = max(abs(cx(i,j,k)), abs(cy(i,j,k))) + 1.-sin_sg(i,j,5)
                  cmax(k) = max( cmax_t, cmax(k) )
               enddo
            enddo
         endif
      enddo
      call mp_reduce_max(cmax,npz)

! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
      c_global = cmax(1)
      if ( npz /= 1 ) then                ! if NOT shallow water test case
         do k=2,npz
            c_global = max(cmax(k), c_global)
         enddo
      endif
      nsplt = int(1. + c_global)
      if ( is_master() .and. nsplt > 3 )  write(*,*) 'Tracer_2d_split=', nsplt, c_global
   else
      nsplt = q_split
   endif
!--------------------------------------------------------------------------------

   frac  = 1. / real(nsplt)

      if( nsplt /= 1 ) then
!$omp parallel do default(shared)
          do k=1,npz
             do j=jsd,jed
                do i=is,ie+1
                   cx(i,j,k) =  cx(i,j,k) * frac
                   xfx(i,j,k) = xfx(i,j,k) * frac
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   mfx(i,j,k) = mfx(i,j,k) * frac
                enddo
             enddo

             do j=js,je+1
                do i=isd,ied
                   cy(i,j,k) =  cy(i,j,k) * frac
                  yfx(i,j,k) = yfx(i,j,k) * frac
                enddo
             enddo

             do j=js,je+1
                do i=is,ie
                  mfy(i,j,k) = mfy(i,j,k) * frac
                enddo
             enddo
          enddo
      endif

!$omp parallel do default(shared)
      do k=1,npz
         do j=js,je
            do i=is,ie
               dp2(i,j,k) = dp1(i,j,k) + (mfx(i,j,k) - mfx(i+1,j,k) +  &
                            mfy(i,j,k) - mfy(i,j+1,k)) * rarea(i,j)
            enddo
         enddo
      enddo

    do it=1,nsplt
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
      call mpp_complete_update_domains(q_pack, q, domain)
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')

!$omp parallel do default(shared) private(ra_x, ra_y, fx, fy)
      do k=1,npz

         do j=jsd,jed
            do i=is,ie
               ra_x(i,j) = area(i,j) + xfx(i,j,k) - xfx(i+1,j,k)
            enddo
         enddo
         do j=js,je
            do i=isd,ied
               ra_y(i,j) = area(i,j) + yfx(i,j,k) - yfx(i,j+1,k)
            enddo
         enddo

         do iq=1,nq
            call fv_tp_2d(q(isd,jsd,k,iq), cx(is,jsd,k), cy(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, mfx=mfx(is,js,k), mfy=mfy(is,js,k))
            do j=js,je
               do i=is,ie
                  q(i,j,k,iq) = ( q(i,j,k,iq)*dp1(i,j,k) + &
                                (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j) )/dp2(i,j,k)
               enddo
            enddo

         enddo
      enddo ! npz

      if ( it /= nsplt ) then
                      call timing_on('COMM_TOTAL')
                          call timing_on('COMM_TRACER')
           q_pack = mpp_start_update_domains(q, domain)
                          call timing_off('COMM_TRACER')
                      call timing_off('COMM_TOTAL')
!$omp parallel do default(shared)
           do k=1,npz
              do j=js,je
                 do i=is,ie
                    dp1(i,j,k) = dp2(i,j,k)
                    dp2(i,j,k) = dp1(i,j,k) + (mfx(i,j,k) - mfx(i+1,j,k) +  &
                                 mfy(i,j,k) - mfy(i,j+1,k))*rarea(i,j)
                 enddo
              enddo
           enddo 
      endif

   enddo  ! nsplt

   if ( id_divg > 0 ) then
        rdt = 1./(frac*dt)

!$omp parallel do default(shared)
        do k=1,npz
        do j=js,je
           do i=is,ie
              dp1(i,j,k) = (xfx(i+1,j,k)-xfx(i,j,k) + yfx(i,j+1,k)-yfx(i,j,k))*rarea(i,j)*rdt
           enddo
        enddo
        enddo
   endif

end subroutine tracer_2d

subroutine tracer_2d_nested(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz,   &
                     nq,  hord, q_split, dt, id_divg, q_pack, z_tracer, k_split, neststruct, parent_grid)

      type(fv_grid_bounds_type), intent(IN) :: bd
      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    ! number of tracers to be advected
      integer, intent(IN) :: hord
      integer, intent(IN) :: q_split, k_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt
      logical, intent(IN) :: z_tracer
      integer, intent(inout) :: q_pack
      real   , intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)   ! Tracers
      real   , intent(INOUT) :: dp1(bd%is:bd%ie,bd%js:bd%je,npz)        ! DELP before dyn_core
      real   , intent(INOUT) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,  npz)    ! Mass Flux X-Dir
      real   , intent(INOUT) :: mfy(bd%is:bd%ie  ,bd%js:bd%je+1,npz)    ! Mass Flux Y-Dir
      real   , intent(INOUT) ::  cx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)  ! Courant Number X-Dir
      real   , intent(INOUT) ::  cy(bd%isd:bd%ied,bd%js :bd%je +1,npz)  ! Courant Number Y-Dir
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(fv_nest_type), intent(INOUT) :: neststruct
      type(fv_atmos_type), intent(INOUT) :: parent_grid
      type(domain2d), intent(INOUT) :: domain

! Local Arrays
      real :: dp2(bd%is:bd%ie,bd%js:bd%je,npz)
#ifdef FLUXBCS
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je  ,npz,nq)
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1,npz,nq)
#else
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je  )
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1)
#endif
      real :: ra_x(bd%is:bd%ie,bd%jsd:bd%jed)
      real :: ra_y(bd%isd:bd%ied,bd%js:bd%je)
      real :: xfx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)
      real :: yfx(bd%isd:bd%ied,bd%js: bd%je+1, npz)
      real :: cmax(npz)
      real :: cmax_t
      real :: c_global
      real :: frac, rdt
      real, parameter :: esl = 1.E-24
      integer :: nsplt, nsplt_parent, msg_split_steps = 1
      integer :: i,j,k,it,iq,n

      real, pointer, dimension(:,:) :: area, rarea
      real, pointer, dimension(:,:,:) :: sin_sg
      real, pointer, dimension(:,:) :: dxa, dya, dx, dy

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

       area => gridstruct%area
      rarea => gridstruct%rarea

      sin_sg => gridstruct%sin_sg
      dxa    => gridstruct%dxa 
      dya    => gridstruct%dya 
      dx     => gridstruct%dx  
      dy     => gridstruct%dy  

!$omp parallel do default(shared)
      do k=1,npz
         do j=jsd,jed
            do i=is,ie+1
               if (cx(i,j,k) > 0.) then
                  xfx(i,j,k) = cx(i,j,k)*dxa(i-1,j)*dy(i,j)*sin_sg(i-1,j,3)
               else
                  xfx(i,j,k) = cx(i,j,k)*dxa(i,j)*dy(i,j)*sin_sg(i,j,1)
               endif
            enddo
         enddo
         do j=js,je+1
            do i=isd,ied
               if (cy(i,j,k) > 0.) then
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j-1)*dx(i,j)*sin_sg(i,j-1,4)
               else
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j)*dx(i,j)*sin_sg(i,j,2)
               endif
            enddo
         enddo
      enddo

!--------------------------------------------------------------------------------
  if ( q_split == 0 ) then
! Determine nsplt

!$omp parallel do default(shared) private(cmax_t )
      do k=1,npz
         cmax(k) = 0.
         if ( k < 4 ) then
! Top layers: C < max( abs(c_x), abs(c_y) )
            do j=js,je
               do i=is,ie
                  cmax_t  = max( abs(cx(i,j,k)), abs(cy(i,j,k)) )
                  cmax(k) = max( cmax_t, cmax(k) )
               enddo
            enddo
         else
            do j=js,je
               do i=is,ie
                  cmax_t  = max(abs(cx(i,j,k)), abs(cy(i,j,k))) + 1.-sin_sg(i,j,5)
                  cmax(k) = max( cmax_t, cmax(k) )
               enddo
            enddo
         endif
      enddo
      call mp_reduce_max(cmax,npz)

! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
      c_global = cmax(1)
      if ( npz /= 1 ) then                ! if NOT shallow water test case
         do k=2,npz
            c_global = max(cmax(k), c_global)
         enddo
      endif
      nsplt = int(1. + c_global)

#ifdef FLUXBCS
      !!*****NOTE*****
      !! If setting the FLUXBCS directive do note that the current
      !!  version of the code requires that Atm, tile, and pelist be
      !!  brought in through use statements. This code will need re
      !! -writing to avoid this.

      !If using flux BCs, nested grid nsplit must be an
      !even multiple of that on the coarse grid
      if (neststruct%do_flux_BCs .or. (gridstruct%nested .and. neststruct%nestbctype > 1)) then
         !Receive from all parent grids
         if (gridstruct%nested) then
            !!NOTE about mpp_recv/mpp_send and scalars:
            !! When passing a scalar, the second argument is not SIZE (which is known to be 1) but a process ID
            call mpp_recv(nsplt_parent,parent_grid%pelist(1))
            nsplt = ceiling(real(nsplt)/real(nsplt_parent))*nsplt
            nsplt = max(nsplt,nsplt_parent)
            msg_split_steps = nsplt/nsplt_parent
         endif
      endif
#endif
      !if ( master )  write(*,*) 'Tracer_2d_split=', nsplt, c_global 
      if ( is_master() .and. nsplt > 3 )  write(*,*) 'Tracer_2d_split=', nsplt, c_global 
   else
      nsplt = q_split
      if (gridstruct%nested .and. neststruct%nestbctype > 1) msg_split_steps = max(q_split/parent_grid%flagstruct%q_split,1)
   endif

#ifdef FLUXBCS
   !Make sure to send to any nested grids which might be expecting a coarse-grid nsplit.
   !(This is outside the if statement since it could be that the coarse grid uses
   !q_split > 0 but the nested grid has q_split = 0)
   if (neststruct%do_flux_BCs .or. (gridstruct%nested .and. neststruct%nestbctype > 1)) then
      if (ANY(neststruct%child_grids) .and. is_master()) then
         do n=1,size(neststruct%child_grids)
            if (neststruct%child_grids(n) .and. Atm(n)%flagstruct%q_split == 0) then
               do i=1,Atm(n)%npes_this_grid
                  call mpp_send(nsplt,Atm(n)%pelist(i))
               enddo
            endif
         enddo
      endif
   endif
#endif FLUXBCS

!--------------------------------------------------------------------------------

   frac  = 1. / real(nsplt)

      if( nsplt /= 1 ) then
!$omp parallel do default(shared)
          do k=1,npz
             do j=jsd,jed
                do i=is,ie+1
                   cx(i,j,k) =  cx(i,j,k) * frac
                   xfx(i,j,k) = xfx(i,j,k) * frac
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   mfx(i,j,k) = mfx(i,j,k) * frac
                enddo
             enddo

             do j=js,je+1
                do i=isd,ied
                   cy(i,j,k) =  cy(i,j,k) * frac
                  yfx(i,j,k) = yfx(i,j,k) * frac
                enddo
             enddo

             do j=js,je+1
                do i=is,ie
                  mfy(i,j,k) = mfy(i,j,k) * frac
                enddo
             enddo
          enddo
      endif

!$omp parallel do default(shared)
      do k=1,npz
         do j=js,je
            do i=is,ie
               dp2(i,j,k) = dp1(i,j,k) + (mfx(i,j,k) - mfx(i+1,j,k) +  &
                            mfy(i,j,k) - mfy(i,j+1,k)) * rarea(i,j)
            enddo
         enddo
      enddo

    do it=1,nsplt
       if ( gridstruct%nested ) then
          neststruct%tracer_nest_timestep = neststruct%tracer_nest_timestep + 1
       end if
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
      call mpp_complete_update_domains(q_pack, q, domain)
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')
	    
      if (gridstruct%nested) then
            do iq=1,nq
                 call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,:,iq), &
                      !0, 0, npx, npy, npz, real(tracer_nest_timestep), real(nsplt), &
                      0, 0, npx, npy, npz, bd, &
                      real(neststruct%tracer_nest_timestep)+real(nsplt*k_split), real(nsplt*k_split), &
                      var_east_t0=neststruct%q_BC(iq)%east_t0, &
                      var_west_t0=neststruct%q_BC(iq)%west_t0, &
                      var_north_t0=neststruct%q_BC(iq)%north_t0, &
                      var_south_t0=neststruct%q_BC(iq)%south_t0, &
                      var_east_t1=neststruct%q_BC(iq)%east_t1, &
                      var_west_t1=neststruct%q_BC(iq)%west_t1, &
                      var_north_t1=neststruct%q_BC(iq)%north_t1, &
                      var_south_t1=neststruct%q_BC(iq)%south_t1, &
                      bctype=neststruct%nestbctype, &
                      nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )
           enddo
      endif


#ifdef FLUXBCS

! NULL !$omp parallel do default(shared) private(ra_x, ra_y, fx, fy)
      do k=1,npz

         do j=jsd,jed
            do i=is,ie
               ra_x(i,j) = area(i,j) + xfx(i,j,k) - xfx(i+1,j,k)
            enddo
         enddo
         do j=js,je
            do i=isd,ied
               ra_y(i,j) = area(i,j) + yfx(i,j,k) - yfx(i,j+1,k)
            enddo
         enddo

         do iq=1,nq
            call fv_tp_2d(q(isd,jsd,k,iq), cx(is,jsd,k), cy(isd,js,k), &
                          npx, npy, hord, fx(is,js,k,iq), fy(is,js,k,iq), xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, mfx=mfx(is,js,k), mfy=mfy(is,js,k))
         enddo
      enddo

      if (neststruct%do_flux_BCs .or. (gridstruct%nested .and. neststruct%nestbctype > 1) ) then

         !call FCT_PD(q,fx,fy,dp1,npx,npy,npz,nq,area, domain)

         call flux_BCs(fx, fy, it, msg_split_steps, npx, npy, npz, nq, q, dp1, dp2, cx, cy, gridstruct%nested, neststruct, parent_grid)
         !call flux_BCs(fx, fy, it, nsplt, npx, npy, npz, nq, q, dp1, dp2, cx, cy)

      endif
!NULL !$omp parallel do default(shared) private(ra_x, ra_y, fx, fy)
      do k=1,npz
         do iq=1,nq

            do j=js,je
               do i=is,ie
                  q(i,j,k,iq) = ( q(i,j,k,iq)*dp1(i,j,k) + &
                                (fx(i,j,k,iq)-fx(i+1,j,k,iq)+fy(i,j,k,iq)-fy(i,j+1,k,iq))*rarea(i,j) )/dp2(i,j,k)
               enddo
            enddo

         enddo
      enddo ! npz

#else
!$omp parallel do default(shared) private(ra_x, ra_y, fx, fy)
      do k=1,npz

         do j=jsd,jed
            do i=is,ie
               ra_x(i,j) = area(i,j) + xfx(i,j,k) - xfx(i+1,j,k)
            enddo
         enddo
         do j=js,je
            do i=isd,ied
               ra_y(i,j) = area(i,j) + yfx(i,j,k) - yfx(i,j+1,k)
            enddo
         enddo

         do iq=1,nq
            call fv_tp_2d(q(isd,jsd,k,iq), cx(is,jsd,k), cy(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, mfx=mfx(is,js,k), mfy=mfy(is,js,k))
            do j=js,je
               do i=is,ie
                  q(i,j,k,iq) = ( q(i,j,k,iq)*dp1(i,j,k) + &
                                (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j) )/dp2(i,j,k)
               enddo
            enddo

         enddo
      enddo ! npz

#endif
      if ( it /= nsplt ) then
                      call timing_on('COMM_TOTAL')
                          call timing_on('COMM_TRACER')
           q_pack = mpp_start_update_domains(q, domain)
                          call timing_off('COMM_TRACER')
                      call timing_off('COMM_TOTAL')
!$omp parallel do default(shared)
           do k=1,npz
              do j=js,je
                 do i=is,ie
                    dp1(i,j,k) = dp2(i,j,k)
                    dp2(i,j,k) = dp1(i,j,k) + (mfx(i,j,k) - mfx(i+1,j,k) +  &
                                 mfy(i,j,k) - mfy(i,j+1,k))*rarea(i,j)
                 enddo
              enddo
           enddo 
      endif
           !Apply nested-grid BCs
           if ( gridstruct%nested ) then
              do iq=1,nq


                 call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,:,iq), &
                      0, 0, npx, npy, npz, bd, &
                      real(neststruct%tracer_nest_timestep), real(nsplt*k_split), &
                      !0, 0, npx, npy, npz, real(tracer_nest_timestep)+real(nsplt), real(nsplt), &
                      var_east_t0=neststruct%q_BC(iq)%east_t0, &
                      var_west_t0=neststruct%q_BC(iq)%west_t0, &
                      var_north_t0=neststruct%q_BC(iq)%north_t0, &
                      var_south_t0=neststruct%q_BC(iq)%south_t0, &
                      var_east_t1=neststruct%q_BC(iq)%east_t1, &
                      var_west_t1=neststruct%q_BC(iq)%west_t1, &
                      var_north_t1=neststruct%q_BC(iq)%north_t1, &
                      var_south_t1=neststruct%q_BC(iq)%south_t1, &
                      bctype=neststruct%nestbctype, &
                      nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )

              end do
           end if


   enddo  ! nsplt

   if ( id_divg > 0 ) then
        rdt = 1./(frac*dt)

!$omp parallel do default(shared)
        do k=1,npz
        do j=js,je
           do i=is,ie
              dp1(i,j,k) = (xfx(i+1,j,k)-xfx(i,j,k) + yfx(i,j+1,k)-yfx(i,j,k))*rarea(i,j)*rdt
           enddo
        enddo
        enddo
   endif

 end subroutine tracer_2d_nested

#ifdef FLUXBCS

!! nestbctype == 2: use just coarse-grid fluxes on nested grid.
!! nestbctype == 3: use incoming fluxes: when flow is into
!!                  nested grid, nested grid uses coarse-grid fluxes,
!!                  and when flow is out of nested grid, coarse grid
!!                  uses nested-grid fluxes.
!! nestbctype == 4: use just nested-grid fluxes on coarse grid

!! Note that to ensure conservation when using the mass-conserving
!! remap BC, the vertical remapping must not change the column's
!! tracer mass. See fv_dynamics.F90 .

!! The new way (to be implemented):
!! 1. on first split timestep, send coarse grid fluxes to nested grid.
!! 2. Split up coarse grid fluxes onto nested-grid faces.
!! 3. Run tracer advection on both grids, accumulating the boundary fluxes
!!    on both grids, for all split timesteps except the last on the coarse grid.
!! 4. If two-way BCs, send accumulated nested-grid fluxes to coarse grid
!! 5. Do correction step on coarse grid: add difference in to fluxes taken on
!!    last coarse-grid split timestep
!!
!! For two-way BCs we do not replace an incoming flux with an outgoing flux.
 !! (Should we also not INCREASE any outgoing flux?)

subroutine flux_BCs(fx, fy, it, msg_split_steps, npx, npy, npz, nq, &
     q, dp1, dp2, cx, cy, nested, neststruct, parent_grid)

  integer, intent(IN) :: it, msg_split_steps, npx, npy, npz, nq

  real, intent(INOUT) :: fx(is:ie+1,js:je  ,npz,nq)
  real, intent(INOUT) :: fy(is:ie , js:je+1,npz,nq)
  real, intent(INOUT) :: q(isd:ied,jsd:jed,npz,nq)
  real, intent(IN) :: dp1(is:ie,js:je,npz)
  real, intent(IN) :: dp2(is:ie,js:je,npz)
  real, intent(IN) :: cx(is:ie+1,jsd:jed  ,npz)  ! Courant Number X-Dir
  real, intent(IN) ::  cy(isd:ied,js :je +1,npz)  ! Courant Number Y-Dir

  logical, intent(IN) :: nested

  type(fv_nest_type), intent(INOUT), target :: neststruct
  type(fv_atmos_type), intent(INOUT) :: parent_grid

  real, save, allocatable, dimension(:,:,:) :: coarse_fluxes_west, coarse_fluxes_east, coarse_fluxes_south, coarse_fluxes_north
  real, allocatable, dimension(:,:,:) :: nested_data_west, nested_data_east, nested_data_south, nested_data_north

  real :: ebuffer(js:je,npz,nq), wbuffer(js:je,npz,nq)
  real :: sbuffer(is:ie,npz,nq), nbuffer(is:ie,npz,nq)

  !Divided fluxes for nested-grid boundary
  real :: fxWd(npy-1,npz,nq)
  real :: fxEd(npy-1,npz,nq)
  real :: fySd(npx-1,npz,nq)
  real :: fyNd(npx-1,npz,nq)
!!$  real :: fxWd(jsd:jed,npz,nq)
!!$  real :: fxEd(jsd:jed,npz,nq)
!!$  real :: fySd(isd:ied,npz,nq)
!!$  real :: fyNd(isd:ied,npz,nq)

  integer n, iq, i, j, k, ic, jc
  integer js1, je1, is1, ie1
  real fluxsplit, corr

  integer, parameter :: iflux_div = 3

  integer, pointer :: refinement, nestbctype, ioffset, joffset
  logical, pointer :: do_flux_BCs, do_2way_flux_BCs

  logical, dimension(:), pointer :: child_grids

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

  !Note that sends and receipts are only done when mod(it,msg_split_steps) == 1 or 0

  !if (.not. concurrent) call mpp_error(FATAL, "Nested flux BCs can only be used with CONCURRENT nesting.")

  !! Set up local pointers
  refinement          => neststruct%refinement
  do_flux_BCs         => neststruct%do_flux_BCs
  do_2way_flux_BCs    => neststruct%do_2way_flux_BCs
  nestbctype          => neststruct%nestbctype

  ioffset     => neststruct%ioffset     
  joffset     => neststruct%joffset     
  child_grids => neststruct%child_grids 

  if (.not. allocated(nest_fx_west_accum)) then

     if (is == 1) then
        allocate(nest_fx_west_accum(js:je,npz,nq))
     else
        allocate(nest_fx_west_accum(1,1,1))
     endif
     if (ie == npx-1) then
        allocate(nest_fx_east_accum(js:je,npz,nq))
     else
        allocate(nest_fx_east_accum(1,1,1))
     endif

     if (js == 1) then
        allocate(nest_fx_south_accum(is:ie,npz,nq))
     else
        allocate(nest_fx_south_accum(1,1,1))
     endif
     if (je == npy-1) then
        allocate(nest_fx_north_accum(is:ie,npz,nq))
     else
        allocate(nest_fx_north_accum(1,1,1))
     endif

  end if

  !NOTE: changing bounds on coarse_fluxes* to include one extra point at each end
  !      changes answers for nestbctype == 2, but not for 3
  if (nested .and. .not. allocated(coarse_fluxes_west) .and. nestbctype /= 4) then
     allocate(coarse_fluxes_west(-1:npy/refinement+2,npz,nq))
     allocate(coarse_fluxes_east(-1:npy/refinement+2,npz,nq))
     allocate(coarse_fluxes_south(-1:npx/refinement+2,npz,nq))
     allocate(coarse_fluxes_north(-1:npx/refinement+2,npz,nq))
!!$     allocate(coarse_fluxes_west(0:npy/refinement+1,npz,nq))
!!$     allocate(coarse_fluxes_east(0:npy/refinement+1,npz,nq))
!!$     allocate(coarse_fluxes_south(0:npx/refinement+1,npz,nq))
!!$     allocate(coarse_fluxes_north(0:npx/refinement+1,npz,nq))
  endif

!!! Flux BCs: Do transfers, if necessary
  if ((do_flux_BCs .or. nested)) then

     if (nested) then

        if ( (mod(it,msg_split_steps) == 1 .or. msg_split_steps == 1) .and. nestbctype /= 4 ) then
        !if (  it == 1 .and. nestbctype /= 4 ) then
           coarse_fluxes_west = 0.
           coarse_fluxes_east = 0.
           coarse_fluxes_south = 0.
           coarse_fluxes_north = 0.

           !Receive coarse-grid fluxes; save between timesteps

                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
           call mpp_sum(coarse_fluxes_west,  size(coarse_fluxes_west),  (/ parent_grid%pelist, pelist /) )
           call mpp_sum(coarse_fluxes_east,  size(coarse_fluxes_east),  (/ parent_grid%pelist, pelist /) )
           call mpp_sum(coarse_fluxes_south, size(coarse_fluxes_south), (/ parent_grid%pelist, pelist /) )
           call mpp_sum(coarse_fluxes_north, size(coarse_fluxes_north), (/ parent_grid%pelist, pelist /) )
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')

        endif

        is1 = max(is,1)
        ie1 = min(ie,npx-1)
        js1 = max(js,1)
        je1 = min(je,npy-1)

        fluxsplit = 1./real(refinement*msg_split_steps)

        !Divide up fluxes depending on the division method
        select case (iflux_div)
        case (1)
           if (is == 1) then
              do iq=1,nq
              do k=1,npz
              do j=max(js,1),min(je,npy-1)
                 fxWd(j,k,iq) = coarse_fluxes_west((j-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif
           if (ie == npx-1) then
              do iq=1,nq
              do k=1,npz
              do j=max(js,1),min(je,npy-1)
                 fxEd(j,k,iq) = coarse_fluxes_east((j-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif
           if (js == 1) then
              do iq=1,nq
              do k=1,npz
              do i=max(is,1),min(ie,npx-1)
                 fySd(i,k,iq) = coarse_fluxes_south((i-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif
           if (je == npy-1) then
              do iq=1,nq
              do k=1,npz
              do i=max(is,1),min(ie,npx-1)
                 fyNd(i,k,iq) = coarse_fluxes_north((i-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif

        case (2)
           if (is == 1) then
              call PLM_flux_division(coarse_fluxes_west, fxWd, npy/refinement, &
                   npy-1, npz, nq, refinement, msg_split_steps)
           endif
           if (ie == npx-1) then
              call PLM_flux_division(coarse_fluxes_east, fxEd, npy/refinement, &
                   npy-1, npz, nq, refinement, msg_split_steps)
           endif
           if (js == 1) then
              call PLM_flux_division(coarse_fluxes_south, fySd, npx/refinement, &
                   npx-1, npz, nq, refinement, msg_split_steps)
           endif
           if (je == npy-1) then
              call PLM_flux_division(coarse_fluxes_north, fyNd, npx/refinement, &
                   npx-1, npz, nq, refinement, msg_split_steps)
           endif

        case (3)
           if (is == 1) then
              call PPM_flux_division(coarse_fluxes_west, fxWd, npy/refinement, &
                   npy-1, npz, nq, refinement, msg_split_steps)
           endif
           if (ie == npx-1) then
              call PPM_flux_division(coarse_fluxes_east, fxEd, npy/refinement, &
                   npy-1, npz, nq, refinement, msg_split_steps)
           endif
           if (js == 1) then
              call PPM_flux_division(coarse_fluxes_south, fySd, npx/refinement, &
                   npx-1, npz, nq, refinement, msg_split_steps)
           endif
           if (je == npy-1) then
              call PPM_flux_division(coarse_fluxes_north, fyNd, npx/refinement, &
                   npx-1, npz, nq, refinement, msg_split_steps)
           endif

        case DEFAULT
           call mpp_error(FATAL, 'iflux_div value not implemented')
        end select

        !replace nested-grid fluxes as desired; note that grids are aligned, making this easier.

        if (nestbctype == 2) then
           if (is == 1) then
              do iq=1,nq
              do k=1,npz
              do j=js1,je1
                 fx(1,j,k,iq) = fxWd(j,k,iq) !coarse_fluxes_west((j-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif
           if (ie == npx-1) then
              do iq=1,nq
              do k=1,npz
              do j=js1,je1
                 fx(npx,j,k,iq) = fxEd(j,k,iq) !coarse_fluxes_east((j-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif

           if (js == 1) then
              do iq=1,nq
              do k=1,npz
              do i=is1,ie1
                 fy(i,1,k,iq) = fySd(i,k,iq) !coarse_fluxes_south((i-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif
           if (je == npy-1) then
              do iq=1,nq
              do k=1,npz
              do i=is1,ie1
                 fy(i,npy,k,iq) = fyNd(i,k,iq) !coarse_fluxes_north((i-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif 
       else if (nestbctype == 3) then
          ! At the boundary, use the flux directed outward from the
          ! grid it is from. In particular a cell's inward flux
          ! (regardless of whether that cell is on the nested or the
          ! coarse grid) should never be replaced by an outward flux,
          ! lest positivity be violated. If BOTH fluxes are directed
          ! inward from their grid, set the boundary flux to
          ! zero. When both are directed outward, the choice is
          ! arbitrary; use the flux of smaller magnitude. If both are in
          ! the same direction we use the smaller of the fluxes.

           if (is == 1) then
              do iq=1,nq
              do k=1,npz
              do j=js,je
                 if (fxWd(j,k,iq)*fx(1,j,k,iq) <= 0.) then
                    !Two incoming fluxes => set to zero
                    !This could most efficiently be represented by
                    ! fx = sign(min(fxWd,-fx,0),-(fxwd+fx)) 
                    ! but I have not yet tested this
                    if (fxWd(j,k,iq) < 0.) then
                       fx(1,j,k,iq) = 0.
                    else
                       fx(1,j,k,iq) = sign(min(fxWd(j,k,iq),-fx(1,j,k,iq)),-(fxWd(j,k,iq)+fx(1,j,k,iq)))
                    endif
                 else
                    !Go with the flux of smaller magnitude
                    fx(1,j,k,iq) = sign(min(abs(fx(1,j,k,iq)),abs(fxWd(j,k,iq))), fx(1,j,k,iq))
                 endif
!!$                 jc = (j-1)/refinement + 1
!!$                 if ( cx(1,j,k) > 0. ) then
!!$                    fx(1,j,k,iq) = fxWd(j,k,iq)
!!$                 endif
              enddo
              enddo
              enddo
           endif
           if (ie == npx-1) then
              do iq=1,nq
              do k=1,npz
              do j=js,je
                 if (fxEd(j,k,iq)*fx(npx,j,k,iq) <= 0.) then
                    if (fxEd(j,k,iq) > 0.) then
                       fx(npx,j,k,iq) = 0.
                    else
                       fx(npx,j,k,iq) = sign(min(-fxEd(j,k,iq),fx(npx,j,k,iq)),-(fxEd(j,k,iq)+fx(npx,j,k,iq)))
                    endif
                 else
                    !Go with the flux of smaller magnitude
                    fx(npx,j,k,iq) = sign(min(abs(fx(npx,j,k,iq)),abs(fxEd(j,k,iq))), fx(npx,j,k,iq))
                 endif
!!$                 jc = (j-1)/refinement + 1
!!$                 if ( cx(npx,j,k) < 0. ) then
!!$                    fx(npx,j,k,iq) = fxEd(j,k,iq)
!!$                 endif
              enddo
              enddo
              enddo
           endif

           if (js == 1) then
              do iq=1,nq
              do k=1,npz
              do i=is,ie
                 if (fySd(i,k,iq)*fy(i,1,k,iq) <= 0.) then
                    if (fySd(i,k,iq) < 0.) then
                       fy(i,1,k,iq) = 0.
                    else
                       fy(i,1,k,iq) = sign(min(fySd(i,k,iq),-fy(i,1,k,iq)),-(fySd(i,k,iq)+fy(i,1,k,iq)))
                    endif
                 else
                    !Go with the flux of smaller magnitude
                    fy(i,1,k,iq) = sign(min(abs(fy(i,1,k,iq)),abs(fySd(i,k,iq))), fy(i,1,k,iq))
                 endif
!!$                 ic = (i-1)/refinement + 1
!!$                 if ( cy(i,1,k) > 0. ) then
!!$                    fy(i,1,k,iq) = fySd(i,k,iq)
!!$                 endif
              enddo
              enddo
              enddo
           endif
           if (je == npy-1) then
              do iq=1,nq
              do k=1,npz
              do i=is,ie
                 if (fyNd(i,k,iq)*fy(i,npy,k,iq) <= 0.) then
                    if (fyNd(i,k,iq) > 0.) then
                       fy(i,npy,k,iq) = 0.
                    else
                       fy(i,npy,k,iq) = sign(min(-fyNd(i,k,iq),fy(i,npy,k,iq)),-(fyNd(i,k,iq)+fy(i,npy,k,iq)))
                    endif
                 else
                    !Go with the flux of smaller magnitude
                    fy(i,npy,k,iq) = sign(min(abs(fy(i,npy,k,iq)),abs(fyNd(i,k,iq))), fy(i,npy,k,iq))
                 endif
!!$                 ic = (i-1)/refinement + 1
!!$                 if ( cy(i,npy,k) < 0. ) then
!!$                    fy(i,npy,k,iq) = fyNd(i,k,iq)
!!$                 endif
              enddo
              enddo
              enddo
           endif            
        else if (nestbctype == 5) then
           !Same as 2, except we do the same test as before: never replace an incoming flux with an outgoing flux
           !ie. if fx(i,j,k) <0. .or. coarse(i,j,k) > 0. then...
           !(Not fully implemented)
           if (is == 1) then
              do iq=1,nq
              do k=1,npz
                 if (js == 1) then
                    if (fx(1,1,k,iq) < 0. .or. coarse_fluxes_west(1,k,iq) > 0.) &
                         fx(1,1,k,iq) =  coarse_fluxes_west(1,k,iq)*fluxsplit
                    if (fy(1,1,k,iq) < 0. .or. coarse_fluxes_south(1,k,iq) > 0.) &
                         fy(1,1,k,iq) =  coarse_fluxes_south(1,k,iq)*fluxsplit
                 endif
                 if (je == npy-1) then
                    if (fx(1,npy-1,k,iq) < 0. .or. coarse_fluxes_west(npy/refinement,k,iq) > 0.) &
                         fx(1,npy-1,k,iq) =  coarse_fluxes_west(npy/refinement,k,iq)*fluxsplit
                    if (fy(1,npy,k,iq) > 0. .or. coarse_fluxes_north(1,k,iq) < 0.) &
                         fy(1,npy,k,iq) =  coarse_fluxes_north(1,k,iq)*fluxsplit
                 endif
              do j=js1,je1
                 if (fx(1,j,k,iq) < 0. .or. coarse_fluxes_west((j-1)/refinement + 1,k,iq) > 0.) &
                         fx(1,j,k,iq) = coarse_fluxes_west((j-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif
           if (ie == npx-1) then
              do iq=1,nq
              do k=1,npz
                 if (js == 1) then
                    if (fx(npx,1,k,iq) > 0. .or. coarse_fluxes_east(1,k,iq) < 0.) &
                         fx(npx,1,k,iq) =  coarse_fluxes_east(1,k,iq)*fluxsplit
                    if (fy(npx-1,1,k,iq) < 0. .or. coarse_fluxes_south(npx/refinement,k,iq) > 0.) &
                         fy(npx-1,1,k,iq) =  coarse_fluxes_south(npx/refinement,k,iq)*fluxsplit
                 endif
                 if (je == npy-1) then
                    if (fx(npx,npy-1,k,iq) > 0. .or. coarse_fluxes_east(npy/refinement,k,iq) < 0.) &
                         fx(npx,npy-1,k,iq) =  coarse_fluxes_east(npy/refinement,k,iq)*fluxsplit
                    if (fy(npx-1,npy,k,iq) > 0. .or. coarse_fluxes_north(npx/refinement,k,iq) < 0.) &
                         fy(npx-1,npy,k,iq) =  coarse_fluxes_north(npx/refinement,k,iq)*fluxsplit
                 endif
              do j=js1,je1
                 if (fx(1,j,k,iq) > 0. .or. coarse_fluxes_west((j-1)/refinement + 1,k,iq) < 0.) &
                         fx(npx,j,k,iq) = coarse_fluxes_east((j-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif

           if (js == 1) then
              do iq=1,nq
              do k=1,npz
              do i=is1,ie1
                 fy(i,1,k,iq) = coarse_fluxes_south((i-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif
           if (je == npy-1) then
              do iq=1,nq
              do k=1,npz
              do i=is1,ie1
                 fy(i,npy,k,iq) = coarse_fluxes_north((i-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif 
        else
           call mpp_error(FATAL, 'nestbctype used is not supported.')
        endif

     endif

  !For TWO-WAY flux BCs:
  !if (nested .and. (nestbctype == 3 .or. nestbctype == 4)) then
  if ( nested ) then
     
     !Nested-grid: ACCUMULATE boundary fluxes at outflow points

     if (mod(it,msg_split_steps) == 1 .or. msg_split_steps == 1) then
        nest_fx_west_accum  = 0.
        nest_fx_east_accum  = 0.
        nest_fx_south_accum = 0.
        nest_fx_north_accum = 0.
     endif

     if (is == 1)     then
        do iq=1,nq
        do k=1,npz
        do j=js,je
           nest_fx_west_accum(j,k,iq)  = nest_fx_west_accum(j,k,iq)  + fx(1,j,k,iq)
        enddo
        enddo
        enddo
     endif
     if (ie == npx-1)     then
        do iq=1,nq
        do k=1,npz
        do j=js,je
           nest_fx_east_accum(j,k,iq)  = nest_fx_east_accum(j,k,iq)  + fx(npx,j,k,iq)
        enddo
        enddo
        enddo
     endif

     if (js == 1)     then
        do iq=1,nq
        do k=1,npz
        do i=is,ie
           nest_fx_south_accum(i,k,iq)  = nest_fx_south_accum(i,k,iq)  + fy(i,1,k,iq)
        enddo
        enddo
        enddo
     endif
     if (je == npy-1)     then
        do iq=1,nq
        do k=1,npz
        do i=is,ie
           nest_fx_north_accum(i,k,iq)  = nest_fx_north_accum(i,k,iq)  + fy(i,npy,k,iq)
        enddo
        enddo
        enddo
     endif

  endif

     !Send coarse-grid boundary fluxes to nested grid
     if (do_flux_BCs .and.  (mod(it,msg_split_steps) == 1 .or. msg_split_steps == 1) .and. nestbctype /= 4) then
        do n=1,size(child_grids)
           if (child_grids(n) .and. nestbctype > 1 .and. nestbctype /= 4) &
                call send_coarse_fluxes(fx, fy, Atm(n)%npx, Atm(n)%npy, ioffset, joffset, &
                refinement, Atm(n)%pelist, Atm(n)%npes_this_grid, npx, npy, npz, nq, neststruct%parent_tile)
        enddo
     endif

  endif

!  if ((do_2way_flux_BCs .or. (nested .and. (nestbctype == 3 .or. nestbctype == 4))) ) then
  if (do_2way_flux_BCs .or. nested ) then

!!$     !Flux check
!!$     if (nested) then
!!$        sumW = 0.
!!$        sumE = 0.
!!$        sumS = 0.
!!$        sumN = 0.
!!$        if (is == 1)     sumW = sum(nest_fx_west_accum(js:je,:,1))
!!$        if (ie == npx-1) sumE = sum(nest_fx_east_accum(js:je,:,1))
!!$        if (js == 1)     sumS = sum(nest_fx_south_accum(is:ie,:,1))
!!$        if (je == npy-1) sumN = sum(nest_fx_north_accum(is:ie,:,1))
!!$        call mpp_sum(sumW)
!!$        call mpp_sum(sumE)
!!$        call mpp_sum(sumS)
!!$        call mpp_sum(sumW)
!!$        if (master) then
!!$           print*, 'Accumulated nested fluxes:'
!!$           print*, 'WEST: ', sumW
!!$           print*, 'EAST: ', sumE
!!$           print*, 'SOUTH: ', sumS
!!$           print*, 'NORTH: ', sumN
!!$        endif
!!$     endif

     if (do_2way_flux_BCs) then
        do n=1,size(child_grids)
           if (child_grids(n) .and. (nestbctype == 3 .or. nestbctype == 4)) then
              !RECEIVE fluxes; note that this will wait until the child grid is on its last split timestep
              !REPLACE fluxes as necessary
              
              call receive_nested_fluxes(fx,fy,Atm(n)%npx, Atm(n)%npy, ioffset, joffset, &
                   refinement, Atm(n)%pelist, Atm(n)%npes_this_grid, &
                   npz, nq, neststruct%parent_tile, nestbctype, q, dp1)
              
           endif
        enddo

     endif

     if (ANY(child_grids)) then

!!$        !Flux check
!!$        sumW = 0.
!!$        sumE = 0.
!!$        sumS = 0.
!!$        sumN = 0.
!!$        if (is == 1)     sumW = sum(nest_fx_west_accum(js:je,:,1))
!!$        if (ie == npx-1) sumE = sum(nest_fx_east_accum(js:je,:,1))
!!$        if (js == 1)     sumS = sum(nest_fx_south_accum(is:ie,:,1))
!!$        if (je == npy-1) sumN = sum(nest_fx_north_accum(is:ie,:,1))
!!$        call mpp_sum(sumW)
!!$        call mpp_sum(sumE)
!!$        call mpp_sum(sumS)
!!$        call mpp_sum(sumW)
!!$        if (master) then
!!$           print*, 'Ending coarse fluxes:'
!!$           print*, 'WEST: ', sumW
!!$           print*, 'EAST: ', sumE
!!$           print*, 'SOUTH: ', sumS
!!$           print*, 'NORTH: ', sumN
!!$        endif

     endif


     if ( nested .and. (nestbctype == 3 .or. nestbctype == 4)  .and. mod(it,msg_split_steps) == 0) then
        call send_nested_fluxes(ioffset, joffset, refinement, parent_grid%pelist, &
             parent_grid%npes_this_grid, npx, npy, npz, nq)
     endif

  endif

end subroutine flux_BCs

subroutine PLM_flux_division(coarse_flux, out_flux, npts_coarse, npts_out, npz, nq, R, split_steps)

  real, intent(IN) :: coarse_flux(-1:npts_coarse+2, npz, nq) !Will need an extra point in both directions
  real, intent(OUT) :: out_flux(npts_out, npz, nq)
  integer, intent(IN) ::  npts_coarse, npts_out, npz, nq, R, split_steps

  integer :: i, n, k, iq

  real :: slope, B, rR, rsplit

  rR = 1./real(R)
  rsplit = 1./real(split_steps)
  
  !Assume evenly spaced cells
  do iq=1,nq
  do k=1,npz
  do i=1,npts_coarse

     !Need to limit slopes to ensure fluxes do not change sign?
!!$     B = coarse_flux(i,k,iq)
!!$     slope = coarse_flux(i+1,k,iq) - coarse_flux(i-1,k,iq)
     B = coarse_flux(i,k,iq)
     slope = min(max(coarse_flux(i+1,k,iq) - coarse_flux(i-1,k,iq),-2.*abs(B)),2.*abs(B))


     do n=0,R-1
        out_flux((i-1)*R+n+1,k,iq) = rsplit * rR * ( 0.5*slope*rR * real(2*n - R + 1)  + B )
     enddo

  enddo
  enddo
  enddo  

end subroutine PLM_flux_division

subroutine PPM_flux_division(coarse_flux, out_flux, npts_coarse, npts_out, npz, nq, R, split_steps)

  real, intent(IN) :: coarse_flux(-1:npts_coarse+2, npz, nq) !Will need an extra point in both directions
  real, intent(OUT) :: out_flux(npts_out, npz, nq)
  integer, intent(IN) ::  npts_coarse, npts_out, npz, nq, R, split_steps

  integer :: i, n, k, iq

  real :: rR, rsplit, a0, a1, a2, viol, scal, d
  real :: fhat(0:npts_coarse) !Interpolated cell-face values: fhat(i) = \hat{f}_{i+1/2}

!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: p1 =  7./12.     ! 0.58333333
  real, parameter:: p2 = -1./12.
  real, parameter:: r3 =  1./3.

  rR = 1./real(R)
  rsplit = 1./real(split_steps)
  
  !Assume evenly spaced cells
  do iq=1,nq
  do k=1,npz

     do i=0,npts_coarse
        fhat(i) = p1*(coarse_flux(i,k,iq)+coarse_flux(i+1,k,iq)) + p2*(coarse_flux(i-1,k,iq)+coarse_flux(i+2,k,iq))
!!$        !Limiting process (see explanation below)
!!$        if (coarse_flux(i,k,iq)*coarse_flux(i+1,k,iq) < 0.) then
!!$           !If signs differ set to zero
!!$           fhat(i) = 0.
!!$        else if (coarse_flux(i,k,iq) > 0.) then
!!$           fhat(i) = max(fhat(i),0.)
!!$        else
!!$           fhat(i) = min(fhat(i),0.)
!!$        endif
     enddo

     do i=1,npts_coarse

        !Construct integral; using my formulation where 0 = midpoint of interval [-1/2,1/2]
        a0 = coarse_flux(i,k,iq)*1.5 - 0.25*(fhat(i)+fhat(i-1))
        a1 = 2.*(fhat(i)-fhat(i-1))
        a2 = 3.*(fhat(i)+fhat(i-1)) - 6.*coarse_flux(i,k,iq)

        !Limit to avoid sign changes
        !Note that this assumes that the reconstruction is zero at
        ! the cell faces between when fluxes change sign. (Something we took care of earlier) We may want
        ! to relax this restriction so that, instead of the entire
        ! reconstructions being positive definite (for positive
        ! fluxes) that merely the divided fluxes are positive
        ! -definite

        !Here, we just need to check whether the extremum is inside the cell and then scale appropriately


        do n=0,R-1
           d = 0.5 - real(n+1)*rR
           out_flux((i-1)*R+n+1,k,iq) = (a0 + a1*d + a2*d**2.)*rR + 0.5*rR**2.*(a1 + 2.*a2*d) + r3*rR**3.*a2
        enddo
        
     enddo
  enddo
  enddo  

end subroutine PPM_flux_division

!do_FCT
!subroutine FCT_PD(q,fx,fy,dp1,npx,npy,npz,nq,wou)
subroutine FCT_PD(q,fx,fy,dp1,npx,npy,npz,nq, area, domain)

  real, intent(inOUT):: q(isd:ied, jsd:jed,npz,nq)
  real, intent(INOUT) :: fx(is:ie+1,js:je  ,npz,nq)
  real, intent(INOUT) :: fy(is:ie , js:je+1,npz,nq)
  real, intent(IN) :: dp1(is:ie,js:je,npz)
  integer, intent(in) :: npx, npy, npz, nq
  real, intent(IN) :: area(isd:ied,jsd:jed)
  type(domain2d), intent(INOUT) :: domain

  !real, dimension(isd:ied,jsd:jed,npz,nq), INTENT(OUT):: wou 
  real, dimension(isd:ied,jsd:jed,npz,nq) :: wou 
  real, parameter:: esl = 1.E-100

  integer i,j,k, n

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed
  
  wou = 1.

  do n=1,nq
     do k=1,npz  

        do j=js,je
           do i=is,ie
              wou(i,j,k,n) = max(0.,fx(i+1,j,k,n)) - min(0.,fx(i,  j,k,n)) +   &
                             max(0.,fy(i,j+1,k,n)) - min(0.,fy(i,  j,k,n)) + esl
              wou(i,j,k,n) = max(0., q(i,j,k,n)*dp1(i,j,k)) / wou(i,j,k,n)*area(i,j)
           enddo
        enddo
     enddo
  enddo

                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
  call mpp_update_domains(wou,domain, complete=.true.)
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')

  do n=1,nq
     do k=1,npz  
      do j=js,je
         do i=is,ie+1
            if ( fx(i,j,k,n) > 0. ) then
               fx(i,j,k,n) = max(0.,min(1., wou(i-1,j,k,n))) * fx(i,j,k,n)
            else
               fx(i,j,k,n) = max(0.,min(1., wou(i,j,k,n))) * fx(i,j,k,n)
            endif
         enddo
      enddo
      do j=js,je+1
         do i=is,ie
            if ( fy(i,j,k,n) > 0. ) then
               fy(i,j,k,n) = max(0.,min(1., wou(i,j-1,k,n))) * fy(i,j,k,n)
            else
               fy(i,j,k,n) = max(0.,min(1., wou(i,j,k,n))) * fy(i,j,k,n)
            endif
         enddo
      enddo
     enddo
  enddo
end subroutine FCT_PD

subroutine send_coarse_fluxes(fx,fy,npx_n,npy_n,ioffset,joffset,refinement, child_pelist, npes_n, npx, npy, npz, nq, tile_with_nest)
  
  real, intent(IN) :: fx(is:ie+1,js:je  ,npz,nq)
  real, intent(IN) :: fy(is:ie , js:je+1,npz,nq)
  integer, intent(IN) :: npx_n, npy_n, refinement, ioffset, joffset, child_pelist(npes_n), npes_n, npx, npy, npz, nq, tile_with_nest

  real, allocatable, dimension(:,:,:) :: coarse_fluxes_west, coarse_fluxes_east, coarse_fluxes_south, coarse_fluxes_north

  integer :: i, j, k, iq
  integer :: instart, inend, jnstart, jnend

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

  !Get starting and stopping indices for nested grid
  instart = ioffset 
  inend   = ioffset + npx_n/refinement - 1
  jnstart = joffset 
  jnend   = joffset + npy_n/refinement - 1

  allocate(coarse_fluxes_south(instart-2:inend+2,npz,nq))
  allocate(coarse_fluxes_north(instart-2:inend+2,npz,nq))

  allocate(coarse_fluxes_west(jnstart-2:jnend+2,npz,nq))
  allocate(coarse_fluxes_east(jnstart-2:jnend+2,npz,nq))

!!$  allocate(coarse_fluxes_south(instart-1:inend+1,npz,nq))
!!$  allocate(coarse_fluxes_north(instart-1:inend+1,npz,nq))
!!$
!!$  allocate(coarse_fluxes_west(jnstart-1:jnend+1,npz,nq))
!!$  allocate(coarse_fluxes_east(jnstart-1:jnend+1,npz,nq))
!!$
  coarse_fluxes_west  = 0.
  coarse_fluxes_east  = 0.
  coarse_fluxes_south = 0.
  coarse_fluxes_north = 0.

  if (tile == tile_with_nest) then

  if (is <= instart .and. ie >= instart) then
     do iq=1,nq
     do k=1,npz
     do j=max(js,jnstart-2),min(je,jnend+2)
        coarse_fluxes_west(j,k,iq) = fx(instart,j,k,iq)
     enddo
     enddo
     enddo
  endif

  if (is <= inend+1 .and. ie >= inend+1) then
!  if (is <= inend+1 .and. (ie >= inend+1 .or. ie == npx-1)) then !Are we sure we have avoided double counting fluxes?
     !if ie == inend+1 then on the adjacent pe is == inend+2 > inend+1 so no double counting!!
     do iq=1,nq
     do k=1,npz
     do j=max(js,jnstart-2),min(je,jnend+2)
        coarse_fluxes_east(j,k,iq) = fx(inend+1,j,k,iq)
     enddo
     enddo
     enddo
  endif


  if (js <= jnstart .and. je >= jnstart) then
     do iq=1,nq
     do k=1,npz
     do i=max(is,instart-2),min(ie,inend+2)
        coarse_fluxes_south(i,k,iq) = fy(i,jnstart,k,iq)
     enddo
     enddo
     enddo
  endif
  if (js <= jnend+1 .and. je >= jnend+1) then
     do iq=1,nq
     do k=1,npz
     do i=max(is,instart-2),min(ie,inend+2)
        coarse_fluxes_north(i,k,iq) = fy(i,jnend+1,k,iq)
     enddo
     enddo
     enddo
  endif

endif

  !Gather the fluxes
  !Using the simplest solution right now

                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
  call mpp_sum(coarse_fluxes_west,  size(coarse_fluxes_west),  (/ pelist, child_pelist /) )
  call mpp_sum(coarse_fluxes_east,  size(coarse_fluxes_east),  (/ pelist, child_pelist /) )
  call mpp_sum(coarse_fluxes_south, size(coarse_fluxes_south), (/ pelist, child_pelist /) )
  call mpp_sum(coarse_fluxes_north, size(coarse_fluxes_north), (/ pelist, child_pelist /) )
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')


  deallocate( coarse_fluxes_west, coarse_fluxes_east, coarse_fluxes_south, coarse_fluxes_north )


end subroutine send_coarse_fluxes

!Sending averages fluxes before sending them, so that we don't have to send as much data
subroutine send_nested_fluxes(ioffset,joffset,refinement, parent_pelist, npes_n, npx, npy, npz, nq)
  
  integer, intent(IN) :: refinement, parent_pelist(npes_n), npes_n, npx, npy, npz, nq, ioffset, joffset

  integer :: i, j, k, iq, n
  real :: val

  real, allocatable, dimension(:,:,:) :: nested_fluxes_south, nested_fluxes_north, nested_fluxes_west, nested_fluxes_east

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

  allocate(nested_fluxes_south((npx-1)/refinement, npz, nq))
  allocate(nested_fluxes_north((npx-1)/refinement, npz, nq))

  allocate(nested_fluxes_east((npy-1)/refinement, npz, nq))
  allocate(nested_fluxes_west((npy-1)/refinement, npz, nq))

  nested_fluxes_west  = 0.
  nested_fluxes_east  = 0.
  nested_fluxes_south = 0.
  nested_fluxes_north = 0.

  if (is == 1) then
     do iq=1,nq
     do k=1,npz
     do j=js,je
        nested_fluxes_west((j-1)/refinement+1,k,iq) = nested_fluxes_west((j-1)/refinement+1,k,iq) + nest_fx_west_accum(j,k,iq)
     enddo
     enddo
     enddo
  endif
  if (ie == npx-1 ) then
     do iq=1,nq
     do k=1,npz
     do j=js,je
        nested_fluxes_east((j-1)/refinement+1,k,iq) = nested_fluxes_east((j-1)/refinement+1,k,iq) + nest_fx_east_accum(j,k,iq)
     enddo
     enddo
     enddo
  endif

  if (js == 1) then
     do iq=1,nq
     do k=1,npz
     do i=is,ie
        nested_fluxes_south((i-1)/refinement+1,k,iq) = nested_fluxes_south((i-1)/refinement+1,k,iq) + nest_fx_south_accum(i,k,iq)
     enddo
     enddo
     enddo
  endif
  if (je == npy-1) then
     do iq=1,nq
     do k=1,npz
     do i=is,ie
        nested_fluxes_north((i-1)/refinement+1,k,iq) = nested_fluxes_north((i-1)/refinement+1,k,iq) + nest_fx_north_accum(i,k,iq)
     enddo
     enddo
     enddo
  endif

  !Gather the fluxes
  !Using the simplest solution right now
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
  call mpp_sum(nested_fluxes_west,  size(nested_fluxes_west),  (/ parent_pelist, pelist /) )
  call mpp_sum(nested_fluxes_east,  size(nested_fluxes_east),  (/ parent_pelist, pelist /) )
  call mpp_sum(nested_fluxes_south, size(nested_fluxes_south), (/ parent_pelist, pelist /) )
  call mpp_sum(nested_fluxes_north, size(nested_fluxes_north), (/ parent_pelist, pelist /) )
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')

  deallocate( nested_fluxes_south, nested_fluxes_north, nested_fluxes_west, nested_fluxes_east )

end subroutine send_nested_fluxes

subroutine receive_nested_fluxes(fx,fy, npx_n, npy_n, ioffset, joffset, refinement, &
     child_pelist, npes_n, npz, nq, tile_with_nest, nestbctype_n, q, dp1)
  
  real, intent(INOUT) :: fx(is:ie+1,js:je  ,npz,nq)
  real, intent(INOUT) :: fy(is:ie , js:je+1,npz,nq)
  integer, intent(IN) :: npx_n, npy_n, refinement, ioffset, joffset, child_pelist(npes_n), &
       npes_n, npz, nq, tile_with_nest, nestbctype_n
  real, intent(IN) :: q(isd:ied,jsd:jed,npz,nq)
  real, intent(IN) :: dp1(is:ie,js:je,npz)

  integer :: i, j, k, iq
  integer :: instart, inend, jnstart, jnend

  real :: outflux, Ry(js-1:je+1), Rx(is-1:ie+1)

  real :: sumN, sumS, sumE, sumW
  real :: maxN, maxS, maxE, maxW

  real, allocatable, dimension(:,:,:) :: nested_fluxes_south, nested_fluxes_north, nested_fluxes_west, nested_fluxes_east

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

  allocate(nested_fluxes_south((npx_n-1)/refinement, npz, nq))
  allocate(nested_fluxes_north((npx_n-1)/refinement, npz, nq))

  allocate(nested_fluxes_east((npy_n-1)/refinement, npz, nq))
  allocate(nested_fluxes_west((npy_n-1)/refinement, npz, nq))

  nested_fluxes_west  = 0.
  nested_fluxes_east  = 0.
  nested_fluxes_south = 0.
  nested_fluxes_north = 0.

  !Gather the fluxes
  !Using the simplest solution right now
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
  call mpp_sum(nested_fluxes_west,  size(nested_fluxes_west),  (/ pelist, child_pelist /) )
  call mpp_sum(nested_fluxes_east,  size(nested_fluxes_east),  (/ pelist, child_pelist /) )
  call mpp_sum(nested_fluxes_south, size(nested_fluxes_south), (/ pelist, child_pelist /) )
  call mpp_sum(nested_fluxes_north, size(nested_fluxes_north), (/ pelist, child_pelist /) )
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')

  sumN = 0.
  sumS = 0.
  sumE = 0.
  sumW = 0.
  maxN = 0.
  maxS = 0.
  maxE = 0.
  maxW = 0.

  if (tile == tile_with_nest) then

  !Get starting and stopping indices for nested grid
  instart = ioffset
  inend   = ioffset + npx_n/refinement - 1
  jnstart = joffset
  jnend   = joffset + npy_n/refinement - 1

  if (nestbctype_n == 2) then
     
     !I don't recall what this section of code is supposed to do.

      if (is <= instart .and. ie > instart) then
         do iq=1,nq
         do k=1,npz
         do j=max(js,jnstart),min(je,jnend)
               fx(instart,j,k,iq) = fx(instart,j,k,iq) + nested_fluxes_west(j-joffset+1,k,iq)
         enddo
         enddo
         enddo
      endif
      if (is <= inend+1 .and. ie >= inend+1) then
         do iq=1,nq
         do k=1,npz
         do j=max(js,jnstart),min(je,jnend)
               fx(inend+1,j,k,iq) = fx(inend+1,j,k,iq) + nested_fluxes_east(j-joffset+1,k,iq)
         enddo
         enddo
         enddo
      endif


      if (js <= jnstart .and. je > jnstart) then
         do iq=1,nq
         do k=1,npz
         do i=max(is,instart),min(ie,inend)
               fy(i,jnstart,k,iq) = fy(i,jnstart,k,iq) + nested_fluxes_south(i-ioffset+1,k,iq)
         enddo
         enddo
         enddo
      endif
      if (js <= jnend+1 .and. je >= jnend+1) then
         do iq=1,nq
         do k=1,npz
         do i=max(is,instart),min(ie,inend)
               fy(i,jnend+1,k,iq) = fy(i,jnend+1,k,iq) + nested_fluxes_north(i-ioffset+1,k,iq)
         enddo
         enddo
         enddo
      endif

   else

      !This outright replaces the coarse fluxes with nested-grid fluxes
      
    if (is <= instart .and. ie >= instart) then
       do iq=1,nq
       do k=1,npz
          !Ry = 1.
       do j=max(js,jnstart),min(je,jnend)
          !if (iq == 1) maxW = max(maxW,(fx(instart,j,k,iq) - nested_fluxes_west(j-joffset+1,k,iq)))
          fx(instart,j,k,iq) = nested_fluxes_west(j-joffset+1,k,iq)
          !if (iq == 1) sumW = sumW + fx(instart,j,k,iq)
       enddo
       enddo
       enddo
    endif
    if (is <= inend+1 .and. ie >= inend+1) then
       do iq=1,nq
       do k=1,npz
          !Ry = 1.
          i = inend+1
       do j=max(js,jnstart),min(je,jnend)
          !if (iq == 1) maxE = max(maxE,fx(i,j,k,iq) - nested_fluxes_east(j-joffset+1,k,iq))
          fx(i,j,k,iq) = nested_fluxes_east(j-joffset+1,k,iq)
          !if (iq == 1) sumE = sumE + fx(i,j,k,iq)
       enddo
       enddo
       enddo
    endif


    if (js <= jnstart .and. je >= jnstart) then
       do iq=1,nq
       do k=1,npz
          !Rx = 1.
       do i=max(is,instart),min(ie,inend)
          !if (iq == 1) maxS = max(maxS, fy(i,jnstart,k,iq) - nested_fluxes_south(i-ioffset+1,k,iq))
          fy(i,jnstart,k,iq) = nested_fluxes_south(i-ioffset+1,k,iq)
          !if (iq == 1) sumS = sumS + fy(i,jnstart,k,iq)
       enddo
       enddo
       enddo
    endif
    if (js <= jnend+1 .and. je >= jnend+1) then
       do iq=1,nq
       do k=1,npz
          j=jnend+1
       do i=max(is,instart),min(ie,inend)
          !if (iq == 1) maxN = max(maxN, fy(i,j,k,iq) - nested_fluxes_north(i-ioffset+1,k,iq))
          fy(i,j,k,iq) = nested_fluxes_north(i-ioffset+1,k,iq)
          !if (iq == 1) sumN = sumN + fy(i,j,k,iq)
       enddo
       enddo
       enddo
    endif

 endif


    !Flux check
    !if (master) then
!!$    if (max(maxW, maxE, maxS, maxN) > 0.) then
!!$       write(gid+100,*) ''
!!$       if (maxW > 0) write(gid+100,*) 'WEST: ', maxW/maxval(abs(fx(instart,:,:,1)))
!!$       if (maxE > 0) write(gid+100,*) 'EAST: ', maxE/maxval(abs(fx(inend+1,:,:,1)))
!!$       if (maxS > 0) write(gid+100,*) 'SOUTH: ', maxS/maxval(abs(fy(:,jnstart,:,1)))
!!$       if (maxN > 0) write(gid+100,*) 'NORTH: ', maxN/maxval(abs(fy(:,jnend+1,:,1)))
!!$    endif
    !endif

endif

  deallocate( nested_fluxes_south, nested_fluxes_north, nested_fluxes_west, nested_fluxes_east )

end subroutine receive_nested_fluxes

#endif FLUXBCS

end module fv_tracer2d_mod
