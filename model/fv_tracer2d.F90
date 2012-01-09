module fv_tracer2d_mod
   use tp_core_mod,       only: fv_tp_2d
   use fv_grid_tools_mod, only: area, rarea, dxa, dya, dx, dy
   use fv_grid_utils_mod, only: sina_u, sina_v, sin_sg
   use fv_mp_mod,         only: gid, domain, mp_reduce_max,   &
                                ng,isd,ied,jsd,jed,is,js,ie,je
   use mpp_domains_mod,   only: mpp_start_update_domains, mpp_complete_update_domains
   use fv_timing_mod,     only: timing_on, timing_off

implicit none
private

public :: tracer_2d, tracer_2d_1L

!---- version number -----
   character(len=128) :: version = '$Id: fv_tracer2d.F90,v 19.0 2012/01/06 19:57:48 fms Exp $'
   character(len=128) :: tagname = '$Name: siena $'

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
         enddo    ! end iq-loop

         if ( it/=nsplt ) then
              i_pack = mpp_start_update_domains(q, domain)
              do j=js,je
                 do i=is,ie
                    dp1(i,j) = dp2(i,j)
                 enddo
              enddo
         endif
     enddo  ! nsplt

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

end module fv_tracer2d_mod
