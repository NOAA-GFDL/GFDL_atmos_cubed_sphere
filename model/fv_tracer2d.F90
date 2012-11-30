module fv_tracer2d_mod
   use tp_core_mod,       only: fv_tp_2d, copy_corners
   use fv_grid_tools_mod, only: area, rarea, dxa, dya, dx, dy, area
   use fv_grid_utils_mod, only: sina_u, sina_v, sin_sg
   use fv_mp_mod,         only: gid, domain, mp_reduce_max,   &
                                ng,isd,ied,jsd,jed,is,js,ie,je, concurrent, masterproc, mp_gather
   use mpp_domains_mod,   only: mpp_start_update_domains, mpp_complete_update_domains, mpp_update_domains, CGRID_NE
   use fv_timing_mod,     only: timing_on, timing_off
   use boundary_mod,      only: nested_grid_BC_apply, nested_grid_BC_apply_intT
   use fv_current_grid_mod,only:nestbctype, nested, q_east, q_north, q_west, q_south, &
           s_weight, tracer_nest_timestep, nsponge, &
           q_east_t0, q_north_t0, q_west_t0, q_south_t0, k_split
   use fv_current_grid_mod,only:parent_grid, child_grids, master, &
        nest_fx_west_accum, nest_fx_east_accum, nest_fx_south_accum, nest_fx_north_accum, ind_update_h, &
        do_flux_BCs, do_2way_flux_BCs, refinement, pelist, tile, ioffset, joffset
   use fv_arrays_mod,     only: Atm
   use mpp_mod,           only: mpp_error, FATAL, mpp_broadcast, mpp_send, mpp_recv, mpp_sum

implicit none
private

public :: tracer_2d, tracer_2d_nested, tracer_2d_1L

!---- version number -----
   character(len=128) :: version = '$Id: fv_tracer2d.F90,v 17.0.2.3.2.13 2012/05/08 20:49:08 Lucas.Harris Exp $'
   character(len=128) :: tagname = '$Name: siena_201211 $'

contains

!-----------------------------------------------------------------------
! !ROUTINE: Perform 2D horizontal-to-lagrangian transport
!-----------------------------------------------------------------------

subroutine tracer_2d_1L(q, dp0, mfx, mfy, cx, cy, npx, npy, npz, nq, hord,  &
                        q_split, k, q3, dt, id_divg)
      integer, intent(IN) :: npx, npy, npz
      integer, intent(IN) :: k
      integer, intent(IN) :: nq    ! number of tracers to be advected
      integer, intent(IN) :: hord
      integer, intent(IN) :: q_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt
      real   , intent(INOUT) :: q(isd:ied,jsd:jed,nq)       ! 2D Tracers
      real   , intent(INOUT) ::q3(isd:ied,jsd:jed,npz,nq)   ! Tracers
      real   , intent(INOUT) :: dp0(is:ie,js:je)        ! DELP before dyn_core
      real   , intent(IN) :: mfx(is:ie+1,js:je)    ! Mass Flux X-Dir
      real   , intent(IN) :: mfy(is:ie  ,js:je+1)    ! Mass Flux Y-Dir
      real   , intent(IN) ::  cx(is:ie+1,jsd:jed)  ! Courant Number X-Dir
      real   , intent(IN) ::  cy(isd:ied,js :je +1)  ! Courant Number Y-Dir

! Local Arrays
      real :: mfx2(is:ie+1,js:je)
      real :: mfy2(is:ie  ,js:je+1)
      real ::  cx2(is:ie+1,jsd:jed)
      real ::  cy2(isd:ied,js :je +1)

      real :: dp1(is:ie,js:je)
      real :: dp2(is:ie,js:je)
      real :: fx(is:ie+1,js:je )
      real :: fy(is:ie , js:je+1)
      real :: ra_x(is:ie,jsd:jed)
      real :: ra_y(isd:ied,js:je)
      real :: xfx(is:ie+1,jsd:jed)
      real :: yfx(isd:ied,js: je+1)
      real :: cmax
      real :: frac, rdt
      integer :: nsplt
      integer :: i,j,it,iq
      integer :: i_pack

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
         if ( gid == 0 .and. nsplt > 5 )  write(6,*) k, 'Tracer_2d_split=', nsplt, cmax
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
                           xfx, yfx, area, ra_x, ra_y, mfx=mfx2, mfy=mfy2 )
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
           if ( nested ) then

              if (concurrent) then
                 call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,iq), &
                      0, 0, npx, npy, real(tracer_nest_timestep)+real(nsplt*k_split), real(nsplt*k_split), &
                      var_east_t0=q_east_t0(:,:,k,iq), &
                      var_west_t0=q_west_t0(:,:,k,iq), &
                      var_north_t0=q_north_t0(:,:,k,iq), &
                      var_south_t0=q_south_t0(:,:,k,iq), &
                      var_east_t1=q_east(:,:,k,iq), &
                      var_west_t1=q_west(:,:,k,iq), &
                      var_north_t1=q_north(:,:,k,iq), &
                      var_south_t1=q_south(:,:,k,iq), &
                      bctype=nestbctype, &
                      nsponge=nsponge, s_weight=s_weight   )
              else
                 call nested_grid_BC_apply(q(isd:ied,jsd:jed,iq), &
                      0, 0, npx, npy, tracer_nest_timestep, nsplt*k_split, &
                      var_east=q_east(:,:,k,iq), &
                      var_west=q_west(:,:,k,iq), &
                      var_north=q_north(:,:,k,iq), &
                      var_south=q_south(:,:,k,iq), bctype=nestbctype, nsponge=nsponge, s_weight=s_weight     )
              endif
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

     if ( nested .and. k == npz) then
        tracer_nest_timestep = tracer_nest_timestep + 1
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


subroutine tracer_2d(q, dp1, mfx, mfy, cx, cy, npx, npy, npz,   &
                     nq,  hord, q_split, dt, id_divg, q_pack, z_tracer)

      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    ! number of tracers to be advected
      integer, intent(IN) :: hord
      integer, intent(IN) :: q_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt
      logical, intent(IN) :: z_tracer
      integer, intent(inout) :: q_pack
      real   , intent(INOUT) :: q(isd:ied,jsd:jed,npz,nq)   ! Tracers
      real   , intent(INOUT) :: dp1(is:ie,js:je,npz)        ! DELP before dyn_core
      real   , intent(INOUT) :: mfx(is:ie+1,js:je,  npz)    ! Mass Flux X-Dir
      real   , intent(INOUT) :: mfy(is:ie  ,js:je+1,npz)    ! Mass Flux Y-Dir
      real   , intent(INOUT) ::  cx(is:ie+1,jsd:jed  ,npz)  ! Courant Number X-Dir
      real   , intent(INOUT) ::  cy(isd:ied,js :je +1,npz)  ! Courant Number Y-Dir

! Local Arrays
      real :: dp2(is:ie,js:je,npz)
      real :: fx(is:ie+1,js:je )
      real :: fy(is:ie , js:je+1)
      real :: ra_x(is:ie,jsd:jed)
      real :: ra_y(isd:ied,js:je)
      real :: xfx(is:ie+1,jsd:jed  ,npz)
      real :: yfx(isd:ied,js: je+1, npz)
      real :: cmax(npz)
      real :: cmax_t
      real :: c_global
      real :: frac, rdt
      integer :: nsplt
      integer :: i,j,k,it,iq

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
      if ( gid == 0 .and. nsplt > 3 )  write(6,*) 'Tracer_2d_split=', nsplt, c_global
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
                          area, ra_x, ra_y, mfx=mfx(is,js,k), mfy=mfy(is,js,k))
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

subroutine tracer_2d_nested(q, dp1, mfx, mfy, cx, cy, npx, npy, npz,   &
                     nq,  hord, q_split, dt, id_divg, q_pack, z_tracer)

      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    ! number of tracers to be advected
      integer, intent(IN) :: hord
      integer, intent(IN) :: q_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt
      logical, intent(IN) :: z_tracer
      integer, intent(inout) :: q_pack
      real   , intent(INOUT) :: q(isd:ied,jsd:jed,npz,nq)   ! Tracers
      real   , intent(INOUT) :: dp1(is:ie,js:je,npz)        ! DELP before dyn_core
      real   , intent(INOUT) :: mfx(is:ie+1,js:je,  npz)    ! Mass Flux X-Dir
      real   , intent(INOUT) :: mfy(is:ie  ,js:je+1,npz)    ! Mass Flux Y-Dir
      real   , intent(INOUT) ::  cx(is:ie+1,jsd:jed  ,npz)  ! Courant Number X-Dir
      real   , intent(INOUT) ::  cy(isd:ied,js :je +1,npz)  ! Courant Number Y-Dir

! Local Arrays
      real :: dp2(is:ie,js:je,npz)
      real :: fx(is:ie+1,js:je  ,npz,nq)
      real :: fy(is:ie , js:je+1,npz,nq)
      real :: ra_x(is:ie,jsd:jed)
      real :: ra_y(isd:ied,js:je)
      real :: xfx(is:ie+1,jsd:jed  ,npz)
      real :: yfx(isd:ied,js: je+1, npz)
      real :: cmax(npz)
      real :: cmax_t
      real :: c_global
      real :: frac, rdt
      real, parameter :: esl = 1.E-24
      integer :: nsplt, nsplt_parent, msg_split_steps = 1
      integer :: i,j,k,it,iq,n

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

      !If using flux BCs, nested grid nsplit must be an
      !even multiple of that on the coarse grid
      if (do_flux_BCs .or. (nested .and. nestbctype > 1)) then
         !Receive from all parent grids
         if (nested) then
            !!NOTE about mpp_recv/mpp_send and scalars:
            !! When passing a scalar, the second argument is not SIZE (which is known to be 1) but a process ID
            call mpp_recv(nsplt_parent,parent_grid%pelist(1))
            nsplt = ceiling(real(nsplt)/real(nsplt_parent))*nsplt
            msg_split_steps = nsplt/nsplt_parent
         endif
      endif
      !if ( gid == masterproc )  write(6,*) 'Tracer_2d_split=', nsplt, c_global 
      if ( gid == masterproc .and. nsplt > 3 )  write(6,*) 'Tracer_2d_split=', nsplt, c_global 
   else
      nsplt = q_split
      if (nested .and. nestbctype > 1) msg_split_steps = q_split/parent_grid%q_split
   endif

   !Make sure to send to any nested grids which might be expecting a coarse-grid nsplit.
   !(This is outside the if statement since it could be that the coarse grid uses
   !q_split > 0 but the nested grid has q_split = 0)
   if (do_flux_BCs .or. (nested .and. nestbctype > 1)) then
      if (ANY(child_grids) .and. masterproc == gid) then
         do n=1,size(child_grids)
            if (child_grids(n) .and. Atm(n)%q_split == 0) then
               do i=1,Atm(n)%npes_this_grid
                  call mpp_send(nsplt,Atm(n)%pelist(i))
               enddo
            endif
         enddo
      endif
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
       if ( nested ) then
          tracer_nest_timestep = tracer_nest_timestep + 1
       end if
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
      call mpp_complete_update_domains(q_pack, q, domain)
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')
	    
              if (nested .and. concurrent) then
            do iq=1,nq
                 call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,:,iq), &
                      !0, 0, npx, npy, npz, real(tracer_nest_timestep), real(nsplt), &
                      0, 0, npx, npy, npz, real(tracer_nest_timestep)+real(nsplt*k_split), real(nsplt*k_split), &
                      var_east_t0=q_east_t0(:,:,:,iq), &
                      var_west_t0=q_west_t0(:,:,:,iq), &
                      var_north_t0=q_north_t0(:,:,:,iq), &
                      var_south_t0=q_south_t0(:,:,:,iq), &
                      var_east_t1=q_east(:,:,:,iq), &
                      var_west_t1=q_west(:,:,:,iq), &
                      var_north_t1=q_north(:,:,:,iq), &
                      var_south_t1=q_south(:,:,:,iq), &
                      bctype=nestbctype, &
                      nsponge=nsponge, s_weight=s_weight   )
           enddo
              endif

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
                          npx, npy, hord, fx(is,js,k,iq), fy(is,js,k,iq), xfx(is,jsd,k), yfx(isd,js,k), &
                          area, ra_x, ra_y, mfx=mfx(is,js,k), mfy=mfy(is,js,k))
         enddo
      enddo


      if (concurrent .and. (do_flux_BCs .or. (nested .and. nestbctype > 1) )) then

         call FCT_PD(q,fx,fy,dp1,npx,npy,npz,nq)

         call flux_BCs(fx, fy, it, msg_split_steps, npx, npy, npz, nq, q, dp1, dp2, cx, cy)

      endif

!$omp parallel do default(shared) private(ra_x, ra_y, fx, fy)
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
           if ( nested ) then
              do iq=1,nq


              if (concurrent) then
                 call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,:,iq), &
                      0, 0, npx, npy, npz, real(tracer_nest_timestep), real(nsplt*k_split), &
                      !0, 0, npx, npy, npz, real(tracer_nest_timestep)+real(nsplt), real(nsplt), &
                      var_east_t0=q_east_t0(:,:,:,iq), &
                      var_west_t0=q_west_t0(:,:,:,iq), &
                      var_north_t0=q_north_t0(:,:,:,iq), &
                      var_south_t0=q_south_t0(:,:,:,iq), &
                      var_east_t1=q_east(:,:,:,iq), &
                      var_west_t1=q_west(:,:,:,iq), &
                      var_north_t1=q_north(:,:,:,iq), &
                      var_south_t1=q_south(:,:,:,iq), &
                      bctype=nestbctype, &
                      nsponge=nsponge, s_weight=s_weight   )
              else
              call nested_grid_BC_apply(q(isd:ied,jsd:jed,:,iq), &
                   0, 0, npx, npy, 1, tracer_nest_timestep, nsplt*k_split, &
                   var_east=q_east(:,:,:,iq), &
                   var_west=q_west(:,:,:,iq), &
                   var_north=q_north(:,:,:,iq), &
                   var_south=q_south(:,:,:,iq), bctype=nestbctype, &
                   nsponge=nsponge, s_weight=s_weight   )
              endif

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

!! nestbctype == 2: use just coarse-grid fluxes on nested grid.
!! nestbctype == 3: use incoming fluxes: when flow is into
!!                  nested grid, nested grid uses coarse-grid fluxes,
!!                  and when flow is out of nested grid, coarse grid
!!                  uses nested-grid fluxes.
!! nestbctype == 4: use just nested-grid fluxes on coarse grid

!! Note that to ensure conservation when using the mass-conserving
!! remap BC, the vertical remapping must not change the column's
!! tracer mass. See fv_dynamics.F90 .

subroutine flux_BCs(fx, fy, it, msg_split_steps, npx, npy, npz, nq, q, dp1, dp2, cx, cy)

  integer, intent(IN) :: it, msg_split_steps, npx, npy, npz, nq

  real, intent(INOUT) :: fx(is:ie+1,js:je  ,npz,nq)
  real, intent(INOUT) :: fy(is:ie , js:je+1,npz,nq)
  real, intent(INOUT) :: q(isd:ied,jsd:jed,npz,nq)
  real, intent(IN) :: dp1(is:ie,js:je,npz)
  real, intent(IN) :: dp2(is:ie,js:je,npz)
  real, intent(IN) :: cx(is:ie+1,jsd:jed  ,npz)  ! Courant Number X-Dir
  real, intent(IN) ::  cy(isd:ied,js :je +1,npz)  ! Courant Number Y-Dir

  real, save, allocatable, dimension(:,:,:) :: coarse_fluxes_west, coarse_fluxes_east, coarse_fluxes_south, coarse_fluxes_north
  real, allocatable, dimension(:,:,:) :: nested_data_west, nested_data_east, nested_data_south, nested_data_north

  real :: ebuffer(js:je,npz,nq), wbuffer(js:je,npz,nq)
  real :: sbuffer(is:ie,npz,nq), nbuffer(is:ie,npz,nq)

  real :: fxd(isd:ied+1,jsd:jed  ,npz,nq)
  real :: fyd(isd:ied  ,jsd:jed+1,npz,nq)

  integer n, iq, i, j, k, ic, jc
  integer js1, je1, is1, ie1
  real fluxsplit, corr

  !Note that sends and receipts are only done when mod(it,msg_split_steps) == 1 or 0

  if (.not. concurrent) call mpp_error(FATAL, "Nested flux BCs can only be used with CONCURRENT nesting.")

  if (nested .and. .not. allocated(coarse_fluxes_west) .and. nestbctype /= 4) then
     allocate(coarse_fluxes_west(npy/refinement,npz,nq))
     allocate(coarse_fluxes_east(npy/refinement,npz,nq))
     allocate(coarse_fluxes_south(npx/refinement,npz,nq))
     allocate(coarse_fluxes_north(npx/refinement,npz,nq))
  endif

!!! Flux BCs: Do transfers, if necessary
  if ((do_flux_BCs .or. nested)) then

     if (nested) then

        if ( (mod(it,msg_split_steps) == 1 .or. msg_split_steps == 1) .and. nestbctype /= 4 ) then
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

        is1 = max(is,2)
        ie1 = min(ie,npx-2)
        js1 = max(js,2)
        je1 = min(je,npy-2)

        fluxsplit = 1./real(refinement*msg_split_steps)
        !replace nested-grid fluxes as desired; note that grids are aligned, making this easier.
        !Also perform FCT so as to ensure positivity is not violated by the replacement fluxes
        !Just correcting for positivity for now
        if (nestbctype == 2) then
           if (is == 1) then
              do iq=1,nq
              do k=1,npz
                 if (js == 1) then
                    fx(1,1,k,iq) =  coarse_fluxes_west(1,k,iq)*fluxsplit
                    fy(1,1,k,iq) =  coarse_fluxes_south(1,k,iq)*fluxsplit
                 endif
                 if (je == npy-1) then
                    fx(1,npy-1,k,iq) =  coarse_fluxes_west(npy/refinement,k,iq)*fluxsplit
                    fy(1,npy,k,iq) =  coarse_fluxes_north(1,k,iq)*fluxsplit
                 endif
              do j=js1,je1
                 fx(1,j,k,iq) = coarse_fluxes_west((j-1)/refinement + 1,k,iq)*fluxsplit
              enddo
              enddo
              enddo
           endif
           if (ie == npx-1) then
              do iq=1,nq
              do k=1,npz
                 if (js == 1) then
                    fx(npx,1,k,iq) =  coarse_fluxes_east(1,k,iq)*fluxsplit
                    fy(npx-1,1,k,iq) =  coarse_fluxes_south(npx/refinement,k,iq)*fluxsplit
                 endif
                 if (je == npy-1) then
                    fx(npx,npy-1,k,iq) =  coarse_fluxes_east(npy/refinement,k,iq)*fluxsplit
                    fy(npx-1,npy,k,iq) =  coarse_fluxes_north(npx/refinement,k,iq)*fluxsplit
                 endif
              do j=js1,je1
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
       else if (nestbctype == 3) then
          ! At the boundary, use the flux directed outward from the
          ! grid it is from. In particular a cell's inward flux
          ! (regardless of whether that cell is on the nested or the
          ! coarse grid) should never be replaced by an outward flux,
          ! lest positivity be violated. If BOTH fluxes are directed
          ! inward from their grid, set the boundary flux to
          ! zero. When both are directed inward, the choice is
          ! arbitrary; here we use the coarse grid's flux

           if (is == 1) then
              do iq=1,nq
              do k=1,npz
              do j=js,je
                 jc = (j-1)/refinement + 1
                 if ( coarse_fluxes_west(jc,k,iq) > 0. ) then
                    fx(1,j,k,iq) = coarse_fluxes_west(jc,k,iq)*fluxsplit
                 else
                    fx(1,j,k,iq) = min(fx(1,j,k,iq),0.)
                 endif
              enddo
              enddo
              enddo
           endif
           if (ie == npx-1) then
              do iq=1,nq
              do k=1,npz
              do j=js,je
                 jc = (j-1)/refinement + 1
                 if ( coarse_fluxes_east(jc,k,iq) < 0. ) then
                    fx(npx,j,k,iq) = coarse_fluxes_east(jc,k,iq)*fluxsplit
                 else
                    fx(npx,j,k,iq) = max(fx(npx,j,k,iq),0.)
                 endif
              enddo
              enddo
              enddo
           endif

           if (js == 1) then
              do iq=1,nq
              do k=1,npz
              do i=is,ie
                 ic = (i-1)/refinement + 1
                 if ( coarse_fluxes_south(ic,k,iq) > 0. ) then
                    fy(i,1,k,iq) = coarse_fluxes_south(ic,k,iq)*fluxsplit
                 else
                    fy(i,1,k,iq) = min(fy(i,1,k,iq),0.)
                 endif
              enddo
              enddo
              enddo
           endif
           if (je == npy-1) then
              do iq=1,nq
              do k=1,npz
              do i=is,ie
                 ic = (i-1)/refinement + 1
                 if ( coarse_fluxes_north(ic,k,iq) < 0. ) then
                    fy(i,npy,k,iq) = coarse_fluxes_north(ic,k,iq)*fluxsplit
                 else
                    fy(i,npy,k,iq) = max(fy(i,npy,k,iq),0.)
                 endif
              enddo
              enddo
              enddo
           endif            
        else if (nestbctype == 4) then
           !Do nothing to coarse grid (Berger and Colella flux BCs)
        else
           call mpp_error(FATAL, 'nestbctype used is not supported.')
        endif

     endif

  !For TWO-WAY flux BCs:
  if (nested .and. (nestbctype == 3 .or. nestbctype == 4)) then
     
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
           if (child_grids(n) .and. Atm(n)%nestbctype > 1 .and. Atm(n)%nestbctype /= 4) &
                call send_coarse_fluxes(fx, fy, Atm(n)%npx, Atm(n)%npy, Atm(n)%ioffset, Atm(n)%joffset, &
                Atm(n)%refinement, Atm(n)%pelist, Atm(n)%npes_this_grid, npx, npy, npz, nq, Atm(n)%parent_tile)
        enddo
     endif

  endif

  if ((do_2way_flux_BCs .or. (nested .and. (nestbctype == 3 .or. nestbctype == 4))) ) then

     if (do_2way_flux_BCs) then
        do n=1,size(child_grids)
           if (child_grids(n) .and. (Atm(n)%nestbctype == 3 .or. Atm(n)%nestbctype == 4)) then
              !RECEIVE fluxes; note that this will wait until the child grid is on its last split timestep
              !REPLACE fluxes as necessary

              call receive_nested_fluxes(fx,fy,Atm(n)%npx, Atm(n)%npy, Atm(n)%ioffset, Atm(n)%joffset, &
                   Atm(n)%refinement, Atm(n)%pelist, Atm(n)%npes_this_grid, &
                   npz, nq, Atm(n)%parent_tile, Atm(n)%nestbctype, q, dp1)
              
           endif
        enddo

     endif

     if (nested .and. (nestbctype == 3 .or. nestbctype == 4)  .and. mod(it,msg_split_steps) == 0) then
        call send_nested_fluxes(ioffset, joffset, refinement, parent_grid%pelist, &
             parent_grid%npes_this_grid, npx, npy, npz, nq)
     endif

  endif

end subroutine flux_BCs

!do_FCT
!subroutine FCT_PD(q,fx,fy,dp1,npx,npy,npz,nq,wou)
subroutine FCT_PD(q,fx,fy,dp1,npx,npy,npz,nq)

  real, intent(inOUT):: q(isd:ied, jsd:jed,npz,nq)
  real, intent(INOUT) :: fx(is:ie+1,js:je  ,npz,nq)
  real, intent(INOUT) :: fy(is:ie , js:je+1,npz,nq)
  real, intent(IN) :: dp1(is:ie,js:je,npz)
  integer, intent(in) :: npx, npy, npz, nq

  !real, dimension(isd:ied,jsd:jed,npz,nq), INTENT(OUT):: wou 
  real, dimension(isd:ied,jsd:jed,npz,nq) :: wou 
  real, parameter:: esl = 1.E-100

  integer i,j,k, n
  
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

  !Get starting and stopping indices for nested grid
  instart = ioffset
  inend   = ioffset + npx_n/refinement - 1
  jnstart = joffset
  jnend   = joffset + npy_n/refinement - 1

  allocate(coarse_fluxes_south(instart:inend,npz,nq))
  allocate(coarse_fluxes_north(instart:inend,npz,nq))

  allocate(coarse_fluxes_west(jnstart:jnend,npz,nq))
  allocate(coarse_fluxes_east(jnstart:jnend,npz,nq))

  coarse_fluxes_west  = 0.
  coarse_fluxes_east  = 0.
  coarse_fluxes_south = 0.
  coarse_fluxes_north = 0.

  if (tile == tile_with_nest) then

  if (is <= instart .and. ie > instart) then
     do iq=1,nq
     do k=1,npz
     do j=max(js,jnstart),min(je,jnend)
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
     do j=max(js,jnstart),min(je,jnend)
        coarse_fluxes_east(j,k,iq) = fx(inend+1,j,k,iq)
     enddo
     enddo
     enddo
  endif


  if (js <= jnstart .and. je > jnstart) then
     do iq=1,nq
     do k=1,npz
     do i=max(is,instart),min(ie,inend)
        coarse_fluxes_south(i,k,iq) = fy(i,jnstart,k,iq)
     enddo
     enddo
     enddo
  endif
  if (js <= jnend+1 .and. je >= jnend+1) then
     do iq=1,nq
     do k=1,npz
     do i=max(is,instart),min(ie,inend)
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

  real, allocatable, dimension(:,:,:) :: nested_fluxes_south, nested_fluxes_north, nested_fluxes_west, nested_fluxes_east

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

  if (tile == tile_with_nest) then

  !Get starting and stopping indices for nested grid
  instart = ioffset
  inend   = ioffset + npx_n/refinement - 1
  jnstart = joffset
  jnend   = joffset + npy_n/refinement - 1

  if (nestbctype_n == 2) then

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

    if (is <= instart .and. ie >= instart) then
       do iq=1,nq
       do k=1,npz
          Ry = 1.
       do j=max(js,jnstart),min(je,jnend)
          fx(instart,j,k,iq) = nested_fluxes_west(j-joffset+1,k,iq)
       enddo
       enddo
       enddo
    endif
    if (is <= inend+1 .and. ie >= inend+1) then
       do iq=1,nq
       do k=1,npz
          Ry = 1.
          i = inend+1
       do j=max(js,jnstart),min(je,jnend)
          fx(i,j,k,iq) = nested_fluxes_east(j-joffset+1,k,iq)
       enddo
       enddo
       enddo
    endif


    if (js <= jnstart .and. je >= jnstart) then
       do iq=1,nq
       do k=1,npz
          Rx = 1.
       do i=max(is,instart),min(ie,inend)
          fy(i,jnstart,k,iq) = nested_fluxes_south(i-ioffset+1,k,iq)
       enddo
       enddo
       enddo
    endif
    if (js <= jnend+1 .and. je >= jnend+1) then
       do iq=1,nq
       do k=1,npz
          j=jnend+1
       do i=max(is,instart),min(ie,inend)
          fy(i,j,k,iq) = nested_fluxes_north(i-ioffset+1,k,iq)
       enddo
       enddo
       enddo
    endif

 endif


endif



  deallocate( nested_fluxes_south, nested_fluxes_north, nested_fluxes_west, nested_fluxes_east )

end subroutine receive_nested_fluxes


end module fv_tracer2d_mod
