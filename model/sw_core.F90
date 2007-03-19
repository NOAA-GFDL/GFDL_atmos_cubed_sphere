      module sw_core

      use mp_mod,     only: ng, is,js,ie,je, isd,jsd,ied,jed,  &
                            mp_corner_comm, tile, domain
      use mpp_domains_mod, only: CGRID_NE, mpp_update_domains
      use grid_tools, only: npx=>npx_g,npy=>npy_g, cosa, sina,  &
                            rdxc, rdyc, dx,dy, dxc,dyc, dxa,dya,  &
                            rdxa, rdya, area, rarea, rarea_c, rdx, rdy
      use tpcore,     only: fv_tp_2d, std_ppm
      use grid_utils, only: edge_vect_s,edge_vect_n,edge_vect_w,edge_vect_e,  &
                            sw_corner, se_corner, ne_corner, nw_corner,       &
                            cosa_u, cosa_v, cosa_s, sina_s, sina_u, sina_v,   &
                            rsin_u, rsin_v, rsin_v, rsina, ec1, ec2, ew, es,  &
                            eww, ess, ue, ve, big_number, da_min_c, fC, f0,   &
                            a11, a12, a21, a22, rsin2, Gnomonic_grid
      use timingModule, only : timing_on, timing_off
#ifdef SW_DYNAMICS
      use test_cases, only: test_case
#endif
      implicit none

! The wave-form is more diffusive than the default "5th" order
      real, parameter:: b1 =   0.0375
      real, parameter:: b2 =  -7./30.
      real, parameter:: b3 =  -23./120.
      real, parameter:: b4 =  13./30.
      real, parameter:: b5 = -11./240.
      real, parameter:: r3  = 1./3.

      private
      public :: c_sw, d_sw, d2a2c_vect

      contains

 
   subroutine c_sw(delpc, delp, ptc, pt, u,v, w, uc,vc, ua,va, wc,  &
                   ut, vt, dt2, hydrostatic)

      real, intent(INOUT), dimension(isd:ied,  jsd:jed+1):: u, vc
      real, intent(INOUT), dimension(isd:ied+1,jsd:jed  ):: v, uc
      real, intent(INOUT), dimension(isd:ied, jsd:jed):: delp,  pt,  ua, va, w
      real, intent(OUT  ), dimension(isd:ied, jsd:jed):: delpc, ptc, ut, vt, wc
      real,    intent(IN) :: dt2
      logical, intent(IN) :: hydrostatic
! Local:
      real, dimension(is-1:ie+1,js-1:je+1):: vort, ke
      real, dimension(is-1:ie+2,js-1:je+1):: fx, fx1, fx2
      real, dimension(is-1:ie+1,js-1:je+2):: fy, fy1, fy2
      real :: dt4
      integer :: i,j, is2, ie1
!                             call timing_on('d2a2c')
      call d2a2c_vect(u, v, ua, va, uc, vc, ut, vt)  ! simpler/quicker version
                                                     ! This version also reduces naturally
                                                     ! for lat-lon or x-y coordinate
!     call d2a2c_vect_v0(u, v, ua, va, uc, vc, ut, vt)
!                             call timing_off('d2a2c')
      do j=js-1,je+1
         do i=is-1,ie+2
            ut(i,j) = dt2*ut(i,j)*dy(i,j)*sina_u(i,j)
         enddo
      enddo
      do j=js-1,je+2
         do i=is-1,ie+1
            vt(i,j) = dt2*vt(i,j)*dx(i,j)*sina_v(i,j)
         enddo
      enddo

!----------------
! Transport delp:
!----------------
! Xdir:
      call fill2_4corners(delp, pt, 1)
      if ( hydrostatic ) then
#ifdef SW_DYNAMICS
           do j=js-1,je+1
              do i=is-1,ie+2      
                 if ( ut(i,j) > 0. ) then
                      fx1(i,j) = delp(i-1,j)
                 else
                      fx1(i,j) = delp(i,j)
                 endif
                 fx1(i,j) =  ut(i,j)*fx1(i,j)
              enddo
           enddo
#else
           do j=js-1,je+1
              do i=is-1,ie+2      
                 if ( ut(i,j) > 0. ) then
                      fx1(i,j) = delp(i-1,j)
                       fx(i,j) =   pt(i-1,j)
                 else
                      fx1(i,j) = delp(i,j)
                       fx(i,j) =   pt(i,j)
                 endif
                 fx1(i,j) =  ut(i,j)*fx1(i,j)
                  fx(i,j) = fx1(i,j)* fx(i,j)
              enddo
           enddo
#endif
      else
           call fill_4corners(w, 1)
           do j=js-1,je+1
              do i=is-1,ie+2      
                 if ( ut(i,j) > 0. ) then
                      fx1(i,j) = delp(i-1,j)
                       fx(i,j) =   pt(i-1,j)
                      fx2(i,j) =    w(i-1,j)
                 else
                      fx1(i,j) = delp(i,j)
                       fx(i,j) =   pt(i,j)
                      fx2(i,j) =    w(i,j)
                 endif
                 fx1(i,j) =  ut(i,j)*fx1(i,j)
                  fx(i,j) = fx1(i,j)* fx(i,j)
                 fx2(i,j) = fx1(i,j)*fx2(i,j)
              enddo
           enddo
      endif

! Ydir:
      call fill2_4corners(delp, pt, 2)
      if ( hydrostatic ) then
           do j=js-1,je+2
              do i=is-1,ie+1      
                 if ( vt(i,j) > 0. ) then
                      fy1(i,j) = delp(i,j-1)
                       fy(i,j) =   pt(i,j-1)
                 else
                      fy1(i,j) = delp(i,j)
                       fy(i,j) =   pt(i,j)
                 endif
                 fy1(i,j) =  vt(i,j)*fy1(i,j)
                  fy(i,j) = fy1(i,j)* fy(i,j)
              enddo
           enddo
           do j=js-1,je+1
              do i=is-1,ie+1    
                 delpc(i,j) = delp(i,j) + (fx1(i,j)-fx1(i+1,j)+fy1(i,j)-fy1(i,j+1))*rarea(i,j)
#ifdef SW_DYNAMICS
                   ptc(i,j) = pt(i,j)
#else
                   ptc(i,j) = (pt(i,j)*delp(i,j) +   &
                              (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j))/delpc(i,j)
#endif
              enddo
           enddo
      else
           call fill_4corners(w, 2)
           do j=js-1,je+2
              do i=is-1,ie+1      
                 if ( vt(i,j) > 0. ) then
                      fy1(i,j) = delp(i,j-1)
                       fy(i,j) =   pt(i,j-1)
                      fy2(i,j) =    w(i,j-1)
                 else
                      fy1(i,j) = delp(i,j)
                       fy(i,j) =   pt(i,j)
                      fy2(i,j) =    w(i,j)
                 endif
                 fy1(i,j) =  vt(i,j)*fy1(i,j)
                  fy(i,j) = fy1(i,j)* fy(i,j)
                 fy2(i,j) = fy1(i,j)*fy2(i,j)
              enddo
           enddo
           do j=js-1,je+1
              do i=is-1,ie+1    
                 delpc(i,j) = delp(i,j) + (fx1(i,j)-fx1(i+1,j)+fy1(i,j)-fy1(i,j+1))*rarea(i,j)
                   ptc(i,j) = (pt(i,j)*delp(i,j) +   &
                              (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j))/delpc(i,j)
                    wc(i,j) = (w(i,j)*delp(i,j) + (fx2(i,j)-fx2(i+1,j) +    &
                               fy2(i,j)-fy2(i,j+1))*rarea(i,j))/delpc(i,j)
              enddo
           enddo
      endif

!------------
! Compute KE:
!------------
#ifdef SMOOTH_GRID
      do j=js-1,je+1
         do i=is-1,ie+1
            if ( ua(i,j) > 0. ) then
                 ke(i,j) = uc(i,j)
            else
                 ke(i,j) = uc(i+1,j)
            endif
         enddo
      enddo
      do j=js-1,je+1
         do i=is-1,ie+1
            if ( va(i,j) > 0. ) then
                 vort(i,j) = vc(i,j)
            else
                 vort(i,j) = vc(i,j+1)
            endif
         enddo
      enddo
#else
      do j=js-1,je+1
         do i=is-1,ie+1
            if ( ua(i,j) > 0. ) then
                 if ( i==1 ) then
                    ke(1,j) = uc(1,  j)*sina_u(1,  j)+v(1,  j)*cosa_u(1,  j)
                 elseif ( i==npx ) then
                    ke(i,j) = uc(npx,j)*sina_u(npx,j)-v(npx,j)*cosa_u(npx,j)
                 else
                    ke(i,j) = uc(i,j)
                 endif
            else
                 if ( i==0 ) then
                    ke(0,j) = uc(1,  j)*sina_u(1,  j)-v(1,  j)*cosa_u(1,  j)
                 elseif ( i==(npx-1) ) then
                    ke(i,j) = uc(npx,j)*sina_u(npx,j)+v(npx,j)*cosa_u(npx,j)
                 else
                    ke(i,j) = uc(i+1,j)
                 endif
            endif
         enddo
      enddo
      do j=js-1,je+1
         do i=is-1,ie+1
            if ( va(i,j) > 0. ) then
               if ( j==1 ) then
                  vort(i,1) = vc(i,  1)*sina_v(i,  1)+u(i,  1)*cosa_v(i,  1)
               elseif ( j==npy ) then
                  vort(i,j) = vc(i,npy)*sina_v(i,npy)-u(i,npy)*cosa_v(i,npy)
               else
                  vort(i,j) = vc(i,j)
               endif
            else
               if ( j==0 ) then
                  vort(i,0) = vc(i,  1)*sina_v(i,  1)-u(i,  1)*cosa_v(i,  1)
               elseif ( j==(npy-1) ) then
                  vort(i,j) = vc(i,npy)*sina_v(i,npy)+u(i,npy)*cosa_v(i,npy)
               else
                  vort(i,j) = vc(i,j+1)
               endif
            endif
         enddo
      enddo
#endif

      dt4 = 0.5*dt2
      do j=js-1,je+1
         do i=is-1,ie+1
            ke(i,j) = dt4*(ua(i,j)*ke(i,j) + va(i,j)*vort(i,j)) 
         enddo
      enddo

!------------------------------
! Compute circulation on C grid
!------------------------------
! To consider using true co-variant winds at face edges?
#ifdef TEST_EDGE
      do j=js-1,je+1
         do i=is,ie+1
            fx(i,j) = uc(i,j) * dxc(i,j)
         enddo
      enddo
      do j=js,je+1
         do i=is-1,ie+1
            fy(i,j) = vc(i,j) * dyc(i,j)
         enddo
      enddo
#else
      is2 = max(2,is); ie1 = min(npx-1,ie+1)
      do j=js-1,je+1
         do i=is2,ie1
            fx(i,j) = uc(i,j)*dxc(i,j)
         enddo
         if(  is   ==  1 ) fx(1,  j) = uc(1,  j)*sina_u(1,  j)*dxc(1,  j)
         if( (ie+1)==npx ) fx(npx,j) = uc(npx,j)*sina_u(npx,j)*dxc(npx,j)
      enddo

      do j=js,je+1
         if( j==1 .or. j==npy ) then
           do i=is-1,ie+1
              fy(i,j) = vc(i,j)*sina_v(i,j)*dyc(i,j)
           enddo
         else
           do i=is-1,ie+1
              fy(i,j) = vc(i,j)*dyc(i,j)
           enddo
         endif
      enddo
#endif
      do j=js,je+1
         do i=is,ie+1
            vort(i,j) =  fx(i,j-1) - fx(i,j) - fy(i-1,j) + fy(i,j)
         enddo
      enddo

! Remove the extra term at the corners:
      if ( sw_corner ) vort(1,    1) = vort(1,    1) + fy(0,   1)
      if ( se_corner ) vort(npx  ,1) = vort(npx,  1) - fy(npx, 1)
      if ( ne_corner ) vort(npx,npy) = vort(npx,npy) - fy(npx,npy)
      if ( nw_corner ) vort(1,  npy) = vort(1,  npy) + fy(0,  npy)

!----------------------------
! Compute absolute vorticity
!----------------------------
      do j=js,je+1
         do i=is,ie+1
#ifdef NON_ROTATING
            vort(i,j) = rarea_c(i,j) * vort(i,j)
#else
            vort(i,j) = fC(i,j) + rarea_c(i,j) * vort(i,j)
#endif
         enddo
      enddo

!----------------------------------
! Transport absolute vorticity:
!----------------------------------
      do j=js,je
         do i=is,ie+1
            if ( i==1 .or. i==npx ) then
                 fy1(i,j) = dt2*v(i,j)*sina_u(i,j)
            else
                 fy1(i,j) = dt2*(v(i,j)-uc(i,j)*cosa_u(i,j))/sina_u(i,j)
            endif
!
            if ( fy1(i,j) > 0. ) then
                 fy(i,j) = vort(i,j)
            else
                 fy(i,j) = vort(i,j+1)
            endif
          enddo
      enddo

      do j=js,je+1
         if ( j==1 .or. j==npy ) then
            do i=is,ie
               fx1(i,j) = dt2*u(i,j)*sina_v(i,j)
               if ( fx1(i,j) > 0. ) then
                    fx(i,j) = vort(i,j)
               else
                    fx(i,j) = vort(i+1,j)
               endif
            enddo
         else
            do i=is,ie
               fx1(i,j) = dt2*(u(i,j)-vc(i,j)*cosa_v(i,j))/sina_v(i,j)
               if ( fx1(i,j) > 0. ) then
                    fx(i,j) = vort(i,j)
               else
                    fx(i,j) = vort(i+1,j)
               endif
            enddo
         endif
      enddo

! Update time-centered winds on the C-Grid
      do j=js,je
         do i=is,ie+1
            uc(i,j) = uc(i,j) + fy1(i,j)*fy(i,j) + rdxc(i,j)*(ke(i-1,j)-ke(i,j))
         enddo
      enddo
      do j=js,je+1
         do i=is,ie
            vc(i,j) = vc(i,j) - fx1(i,j)*fx(i,j) + rdyc(i,j)*(ke(i,j-1)-ke(i,j))
         enddo
      enddo

   end subroutine c_sw



!-------------------------------------------------------------------------------
!
!     d_sw :: D-Grid Shallow Water Routine
!
   subroutine d_sw(delpc, delp,  ptc,   pt, u,  v, w, uc,vc, &
                   ua,va, xflux, yflux, cx, cy,              &
                   crx_adv, cry_adv,  xfx_adv, yfx_adv,      &
                   dt, hord_mt, hord_vt, hord_tm,    &
                   dddmp, hydrostatic, uniform_ppm)

      integer, intent(IN):: hord_mt, hord_vt, hord_tm
      real   , intent(IN):: dt, dddmp
      real, intent(INOUT), dimension(isd:ied,  jsd:jed):: delp, pt, ua, va, w
      real, intent(INOUT), dimension(isd:ied  ,jsd:jed+1):: u, vc
      real, intent(INOUT), dimension(isd:ied+1,jsd:jed  ):: v, uc
      real, intent(OUT),   dimension(isd:ied,  jsd:jed)  :: delpc, ptc
      real, intent(INOUT):: xflux(is:ie+1,js:je  )
      real, intent(INOUT):: yflux(is:ie  ,js:je+1)
      real, intent(INOUT)::    cx(is:ie+1,jsd:jed  )
      real, intent(INOUT)::    cy(isd:ied,js:je+1)
      logical, intent(IN):: uniform_ppm
      logical, intent(IN):: hydrostatic
      real, intent(OUT), dimension(is:ie+1,jsd:jed):: crx_adv, xfx_adv
      real, intent(OUT), dimension(isd:ied,js:je+1):: cry_adv, yfx_adv
! Local:
!     real :: ut(isd:ied+1,jsd:jed  )
!     real :: vt(isd:ied  ,jsd:jed+1)
      real :: ut(is-1:ie+2,jsd: jed  )
      real :: vt(isd: ied ,js-1:je+2)

      real, dimension(is:ie+1,js:je+1):: ub, vb, cfl, fy_ke
      real :: vort(isd:ied,jsd:jed)     ! Vorticity
      real ::   ke(isd:ied+1,jsd:jed+1) ! Kinetic Energy
      real ::   fx(is:ie+1,js:je  )  ! 1-D X-direction Fluxes
      real ::   fy(is:ie  ,js:je+1)  ! 1-D Y-direction Fluxes
      real ::   gy(is:ie  ,js:je+1)  ! work Y-dir flux array

      real :: dt5, dt6, rdt
      real :: utmp, vtmp, damp
      real :: u1, v1, vs, uw, vn, ue
      real :: k1, k2, k3
      integer :: i,j, is2, ie1

      rdt = 1./ dt

#ifdef SW_DYNAMICS
      if ( test_case == 1 ) then
        do j=jsd,jed
           do i=is,ie+1
              xfx_adv(i,j) = dt * uc(i,j) / sina_u(i,j)
              if (xfx_adv(i,j) > 0.) then
                  crx_adv(i,j) = xfx_adv(i,j) * rdxa(i-1,j)
              else
                  crx_adv(i,j) = xfx_adv(i,j) * rdxa(i,j)
              endif
              xfx_adv(i,j) = dy(i,j)*xfx_adv(i,j)*sina_u(i,j)
           enddo
        enddo

        do j=js,je+1
           do i=isd,ied
              yfx_adv(i,j) = dt * vc(i,j) / sina_v(i,j)
              if (yfx_adv(i,j) > 0.) then
                 cry_adv(i,j) = yfx_adv(i,j) * rdya(i,j-1)
              else
                 cry_adv(i,j) = yfx_adv(i,j) * rdya(i,j)
              endif
              yfx_adv(i,j) = dx(i,j)*yfx_adv(i,j)*sina_v(i,j)
           enddo
        enddo
      else
#endif

        do j=jsd,jed
           if(j/=0 .and. j/=1 .and. j/=(npy-1) .and. j/=npy) then
           do i=is-1,ie+2
              vtmp = 0.25*(vc(i-1,j) + vc(i,j) + vc(i-1,j+1) + vc(i,j+1))
              ut(i,j) = (uc(i,j) - cosa_u(i,j)*vtmp) * rsin_u(i,j)
           enddo
! Fix i=1 & npx:
             if(     is==1   ) ut(1  ,j) = uc(1,  j) * rsin_u(1,  j)
             if( (ie+1)==npx ) ut(npx,j) = uc(npx,j) * rsin_u(npx,j)
           endif
        enddo

        do j=js-1,je+2
           if( j/=1 .and. j/=npy ) then
              do i=isd,ied
                 utmp = 0.25*(uc(i,j-1) + uc(i+1,j-1) + uc(i,j) + uc(i+1,j))
                 vt(i,j) = (vc(i,j) - cosa_v(i,j)*utmp) * rsin_v(i,j)
              enddo
           endif
        enddo

! Fix the 4 edges (if existed):
        if ( is==1 ) then
           do j=jsd,jed
              ut(1,j) = uc(1,j) * rsin_u(1,j)
           enddo
           do j=max(3,js), min(npy-2,je+1)
!             vt(0,j) = vc(0,j) - 0.25*cosa_v(0,j)*   &
              vt(0,j) = vc(0,j) + 0.25*cosa_v(1,j)*   &
                                 (ut(0,j)+ut(1,j)+ut(0,j-1)+ut(1,j-1))
              vt(1,j) = vc(1,j) - 0.25*cosa_v(1,j)*   &
                                 (ut(1,j)+ut(2,j)+ut(1,j-1)+ut(2,j-1))
           enddo
        endif

        if ( (ie+1)==npx ) then
           do j=jsd,jed
              ut(npx,j) = uc(npx,j) * rsin_u(npx,j)
           enddo
           do j=max(3,js), min(npy-2,je+1)
              vt(npx-1,j) = vc(npx-1,j) - 0.25*cosa_v(npx-1,j)*   &
                                 (ut(npx-1,j)+ut(npx,j)+ut(npx-1,j-1)+ut(npx,j-1))
!             vt(npx,j) = vc(npx,j) - 0.25*cosa_v(npx,j)*   &
              vt(npx,j) = vc(npx,j) + 0.25*cosa_v(npx-1,j)*   &
                                 (ut(npx,j)+ut(npx+1,j)+ut(npx,j-1)+ut(npx+1,j-1))
           enddo
        endif

        if ( js==1 ) then
           do i=isd,ied
              vt(i,1) = vc(i,1) * rsin_v(i,1)
           enddo
           do i=max(3,is),min(npx-2,ie+1)
!             ut(i,0) = uc(i,0) - 0.25*cosa_u(i,0)*   &
              ut(i,0) = uc(i,0) + 0.25*cosa_u(i,1)*   &
                                 (vt(i,0)+vt(i,1)+vt(i-1,0)+vt(i-1,1))
              ut(i,1) = uc(i,1) - 0.25*cosa_u(i,1)*   &
                                 (vt(i,1)+vt(i,2)+vt(i-1,1)+vt(i-1,2))
           enddo
        endif

        if ( (je+1)==npy ) then
           do i=isd,ied
              vt(i,npy) = vc(i,npy) * rsin_v(i,npy)
           enddo
           do i=max(3,is),min(npx-2,ie+1)
              ut(i,npy-1) = uc(i,npy-1) - 0.25*cosa_u(i,npy-1)*   &
                           (vt(i,npy-1)+vt(i,npy)+vt(i-1,npy-1)+vt(i-1,npy))
!             ut(i,npy) = uc(i,npy) - 0.25*cosa_u(i,npy)*   &
              ut(i,npy) = uc(i,npy) + 0.25*cosa_u(i,npy-1)*   &
                         (vt(i,npy)+vt(i,npy+1)+vt(i-1,npy)+vt(i-1,npy+1))
           enddo
        endif

! Finally, fix the 4 corners:
        if( sw_corner ) then
            utmp = uc(2,1)+uc(2,2)
            vtmp = vc(1,2)+vc(2,2)
            u1 = (utmp - vtmp*cosa(2,2))*rsina(2,2)
            v1 = (vtmp - utmp*cosa(2,2))*rsina(2,2)
            vs = vt(1,1)+vt(2,1)
            uw = ut(1,1)+ut(1,2)
            ut(2,1) = uc(2,1) - 0.25*(v1+vs)*cosa_u(2,1)
            vt(1,2) = vc(1,2) - 0.25*(u1+uw)*cosa_v(1,2) 
! South
            v1 = (vc(1,0)+vc(2,0)+(uc(2,-1)+uc(2,0))*cosa(2,2))*rsina(2,2)
            ut(2,0) = uc(2,0) + 0.25*(v1+vs)*cosa_u(2,1)
! West
            u1 = (uc(0,1)+uc(0,2)+(vc(-1,2)+vc(0,2))*cosa(2,2))*rsina(2,2)
            vt(0,2) = vc(0,2) + 0.25*(u1+uw)*cosa_v(1,2) 
        endif

        if( se_corner ) then
            utmp = uc(npx-1,1)+uc(npx-1,2)
            vtmp = vc(npx-2,2)+vc(npx-1,2)
            u1 = (utmp - vtmp*cosa(npx-1,2))*rsina(npx-1,2)
            v1 = (vtmp - utmp*cosa(npx-1,2))*rsina(npx-1,2)
            vs = vt(npx-1,1)+vt(npx-2,1)
            ue = ut(npx,1)+ut(npx,2)
            ut(npx-1,1) = uc(npx-1,1) - 0.25*(v1+vs)*cosa_u(npx-1,1)
            vt(npx-1,2) = vc(npx-1,2) - 0.25*(u1+ue)*cosa_v(npx-1,2) 
! South
            v1 = (vc(npx-2,0)+vc(npx-1,0) +        &
                 (uc(npx-1,-1)+uc(npx-1,0))*cosa(npx-1,2))*rsina(npx-1,2)
            ut(npx-1,0) = uc(npx-1,0) + 0.25*(v1+vs)*cosa_u(npx-1,1)
! East
            u1 = (uc(npx+1,1)+uc(npx+1,2) +     &
                 (vc(npx  ,2)+vc(npx+1,2))*cosa(npx-1,2))*rsina(npx-1,2)
            vt(npx,2) = vc(npx,2) + 0.25*(u1+ue)*cosa_v(npx-1,2) 
        endif

        if( ne_corner ) then
            utmp = uc(npx-1,npy-2)+uc(npx-1,npy-1)
            vtmp = vc(npx-2,npy-1)+vc(npx-1,npy-1)
            u1 = (utmp - vtmp*cosa(npx-1,npy-1))*rsina(npx-1,npy-1)
            v1 = (vtmp - utmp*cosa(npx-1,npy-1))*rsina(npx-1,npy-1)
            vn = vt(npx-1,npy)+vt(npx-2,npy)
            ue = ut(npx,npx-1)+ut(npx,npx-2)
            ut(npx-1,npx-1) = uc(npx-1,npy-1) - 0.25*(v1+vn)*cosa_u(npx-1,npy-1)
            vt(npx-1,npy-1) = vc(npx-1,npy-1) - 0.25*(u1+ue)*cosa_v(npx-1,npy-1)
! east
            u1 = (uc(npx+1,npy-2)+uc(npx+1,npy-1) +    &
             (vc(npx,npy-1)+vc(npx+1,npy-1))*cosa(npx-1,npy-1))*rsina(npx-1,npy-1)
            vt(npx,npy-1) = vc(npx,npy-1) + 0.25*(u1+ue)*cosa_v(npx-1,npy-1)
! North
            v1 = (vc(npx-2,npy+1)+vc(npx-1,npy+1) +    &
              (uc(npx-1,npy  )+uc(npx-1,npy+1))*cosa(npx-1,npy-1))*rsina(npx-1,npy-1)
            ut(npx-1,npy) = uc(npx-1,npy) + 0.25*(v1+vn)*cosa_u(npx-1,npy-1)
        endif

        if( nw_corner ) then
            utmp = uc(2,npy-1)+uc(2,npy-2)
            vtmp = vc(1,npy-1)+vc(2,npy-1)
            u1 = (utmp - vtmp*cosa(2,npy-1))*rsina(2,npy-1)
            v1 = (vtmp - utmp*cosa(2,npy-1))*rsina(2,npy-1)
            vn = vt(1,npy)+vt(2,npy)
            uw = ut(1,npy-1)+ut(1,npy-2)
            ut(2,npy-1) = uc(2,npy-1) - 0.25*(v1+vn)*cosa_u(2,npy-1)
            vt(1,npy-1) = vc(1,npy-1) - 0.25*(u1+uw)*cosa_v(1,npy-1) 
! North:
            v1 = (vc(1,npy+1)+vc(2,npy+1) +   &
                 (uc(2,npy+1)+uc(2,npy))*cosa(2,npy-1))*rsina(2,npy-1)
            ut(2,npy) = uc(2,npy) + 0.25*(v1+vn)*cosa_u(2,npy-1)
! West:
            u1 = (uc( 0,npy-1)+uc(0,npy-2) +    &
                 (vc(-1,npy-1)+vc(0,npy-1))*cosa(2,npy-1))*rsina(2,npy-1)
            vt(0,npy-1) = vc(0,npy-1) + 0.25*(u1+uw)*cosa_v(1,npy-1) 
        endif
 
        do j=jsd,jed
           do i=is,ie+1
              xfx_adv(i,j) = dt*ut(i,j)
           enddo
        enddo

        do j=js,je+1
           do i=isd,ied
              yfx_adv(i,j) = dt*vt(i,j)
           enddo
        enddo

! Compute E-W CFL number:
        do j=jsd,jed
           do i=is,ie+1
              if (xfx_adv(i,j) > 0.) then
                  crx_adv(i,j) = xfx_adv(i,j) * rdxa(i-1,j)
              else
                  crx_adv(i,j) = xfx_adv(i,j) * rdxa(i,j)
              endif
           enddo
        enddo
        do j=jsd,jed
           do i=is,ie+1
              xfx_adv(i,j) = dy(i,j)*xfx_adv(i,j)*sina_u(i,j)
           enddo
        enddo


        do j=js,je+1
           do i=isd,ied
              if (yfx_adv(i,j) > 0.) then
                 cry_adv(i,j) = yfx_adv(i,j) * rdya(i,j-1)
              else
                 cry_adv(i,j) = yfx_adv(i,j) * rdya(i,j)
              endif
           enddo
        enddo
        do j=js,je+1
           do i=isd,ied
              yfx_adv(i,j) = dx(i,j)*yfx_adv(i,j)*sina_v(i,j)
           enddo
        enddo

#ifdef SW_DYNAMICS
      endif
#endif


! Transport delp
      call fv_tp_2d(delp, crx_adv, cry_adv, npx-1, npy-1, hord_tm,  &
                    fx, fy, xfx_adv,yfx_adv, area, uniform_ppm)

#ifdef SW_DYNAMICS
        do j=js,je
           do i=is,ie
              delp(i,j) = delp(i,j) +    &
                         (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j)
              ptc(i,j) = pt(i,j)
           enddo
        enddo
#else

! <<< Save the mass fluxes to the "Flux Capacitor" for tracer transport >>>
        do j=jsd,jed
           do i=is,ie+1
              cx(i,j) = cx(i,j) + crx_adv(i,j)
           enddo
        enddo       
        do j=js,je
           do i=is,ie+1
              xflux(i,j) = xflux(i,j) + fx(i,j)
           enddo
        enddo       

        do j=js,je+1
           do i=isd,ied
              cy(i,j) = cy(i,j) + cry_adv(i,j)
           enddo
           do i=is,ie
              yflux(i,j) = yflux(i,j) + fy(i,j)
           enddo
        enddo 

        if ( .not. hydrostatic ) then
            call fv_tp_2d(w, crx_adv,cry_adv, npx-1, npy-1, hord_vt, &
                      ub, gy, xfx_adv,yfx_adv, area, uniform_ppm, mfx=fx, mfy=fy)
            do j=js,je
               do i=is,ie
                  w(i,j) = w(i,j)*delp(i,j) +             &
                           (ub(i,j)-ub(i+1,j)+gy(i,j)-gy(i,j+1))*rarea(i,j)
               enddo
            enddo
        endif

        call fv_tp_2d(pt, crx_adv,cry_adv, npx-1, npy-1, hord_tm,  &
                      ub, gy, xfx_adv,yfx_adv, area, uniform_ppm, mfx=fx, mfy=fy)
        do j=js,je
           do i=is,ie
              pt(i,j) = pt(i,j)*delp(i,j) +             &
                       (ub(i,j)-ub(i+1,j)+gy(i,j)-gy(i,j+1))*rarea(i,j)
              delp(i,j) = delp(i,j) +                     &
                         (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j)
              pt(i,j) = pt(i,j) / delp(i,j)
           enddo
        enddo

        if ( .not. hydrostatic ) then
            do j=js,je
               do i=is,ie
                  w(i,j) = w(i,j) / delp(i,j)
               enddo
            enddo
        endif
#endif

#ifdef SW_DYNAMICS
      if (test_case > 1) then
#endif

!----------------------
! Kinetic Energy Fluxes
!----------------------
! Compute B grid contra-variant components for KE:

      do j=js,je+1
         do i=is,ie+1
            cfl(i,j) = 0.5*(cry_adv(i-1,j)+cry_adv(i,j))
         enddo
      enddo
      do j=js,je+1
         do i=is,ie+1
            if ( cfl(i,j) > 0. ) then
                 vb(i,j) = cfl(i,j)*dy(i,j-1)*rdt
            else
                 vb(i,j) = cfl(i,j)*dy(i,j)*rdt
            endif
         enddo
      enddo

      call ytp_v(cfl, u, v, fy_ke, vb, hord_mt)

      do j=js,je+1
         do i=is,ie+1
            cfl(i,j) = 0.5*(crx_adv(i,j-1)+crx_adv(i,j))
         enddo
      enddo
      do j=js,je+1
         do i=is,ie+1
            if ( cfl(i,j) > 0. ) then
                 ub(i,j) = cfl(i,j)*dx(i-1,j)*rdt
            else
                 ub(i,j) = cfl(i,j)*dx(i,j)*rdt
            endif
         enddo
      enddo

      call xtp_u(cfl, u, v, vb, ub, hord_mt)

      dt5 = 0.5*dt
      do j=js,je+1
         do i=is,ie+1
            ke(i,j) = max(0., dt5*(vb(i,j) + fy_ke(i,j)))
         enddo
      enddo

!-----------------------------------------
! Fix KE at the 4 corners of the face:
!-----------------------------------------
!  call timing_on('CORNER_KE')
   if ( Gnomonic_grid ) then
      dt6 = dt / 6.
      if ( sw_corner ) then
           k1 = (uc(1,1) + uc(1,0)) * u(1,1)
           k2 = (vc(1,1) + vc(0,1)) * v(1,1)
           k3 = (uc(1,1) + vc(1,1)) * u(0,1)
           ke(1,1) = dt6*max(0.,k1 + k2 + k3)*rsin_v(1,1)
      endif
      if ( se_corner ) then
           i = npx;      j = 1
           k1 = (uc(i,j) + uc(i,j-1)) * u(i-1,j)
           k2 = (vc(i,j) + vc(i-1,j)) * v(i,j)
           k3 = (uc(i,j) - vc(i-1,j)) * u(i,j)
           ke(i,j) = dt6*max(0.,k1 + k2 + k3)*rsin_v(npx-1,1)
      endif
      if ( ne_corner ) then
           i = npx;      j = npy
           k1 = (uc(i,  j) + uc(i,j-1)) * u(i-1,j)
           k2 = (vc(i,  j) + vc(i-1,j)) * v(i,j-1)
           k3 = (uc(i,j-1) + vc(i-1,j)) * u(i,j)
           ke(i,j) = dt6*max(0.,k1 + k2 + k3)*rsin_v(npx-1,npy)
      endif
      if ( nw_corner ) then
           i = 1;      j = npy
           k1 = (uc(i,  j) + uc(i,j-1)) * u(i,j)
           k2 = (vc(i,  j) + vc(i-1,j)) * v(i,j-1)
           k3 = (uc(i,j-1) - vc(i,  j)) * u(i-1,j)
           ke(i,j) = dt6*max(0.,k1 + k2 + k3)*rsin_v(1,npy)
      endif
   else
      call mp_corner_comm(ke, npx, npy) 
      if (sw_corner) ke(1,    1) = r3*(ke(2,      1)+ke(1,      2)+ke(0,      1))
      if (se_corner) ke(npx,  1) = r3*(ke(npx+1,  1)+ke(npx,    2)+ke(npx-1,  1))
      if (ne_corner) ke(npx,npy) = r3*(ke(npx+1,npy)+ke(npx,npy-1)+ke(npx-1,npy))
      if (nw_corner) ke(1,  npy) = r3*(ke(2,    npy)+ke(1,  npy-1)+ke(0,    npy))
   endif
!  call timing_off('CORNER_KE')

!-----------------------------
! Compute divergence damping
!-----------------------------
!                          dddmp must be < 0.25 for stability
!                          area ~ dxb*dyb*sin(alpha)
          do j=js,je+1
             if ( j==1 .or. j==npy ) then
                do i=is-1,ie+1
                   ptc(i,j) = u(i,j)*dyc(i,j)*sina_v(i,j)
                enddo
             else
                do i=is-1,ie+1
                   ptc(i,j) = (u(i,j)-0.5*(va(i,j-1)+va(i,j))*cosa_v(i,j))   &
                            *dyc(i,j)*sina_v(i,j)
                enddo
             endif
          enddo

          is2 = max(2,is); ie1 = min(npx-1,ie+1)
          do j=js-1,je+1
             do i=is2, ie1
                vort(i,j) = (v(i,j) - 0.5*(ua(i-1,j)+ua(i,j))*cosa_u(i,j))  &
                            *dxc(i,j)*sina_u(i,j)
             enddo
             if (  is   ==  1 ) vort(1,  j) = v(1,  j)*dxc(1,  j)*sina_u(1,  j)
             if ( (ie+1)==npx ) vort(npx,j) = v(npx,j)*dxc(npx,j)*sina_u(npx,j)
          enddo

          do j=js,je+1
             do i=is,ie+1
                delpc(i,j) = -vort(i,j-1) + vort(i,j) - ptc(i-1,j) + ptc(i,j)
             enddo
          enddo

! Remove the extra term at the corners:
          if (sw_corner) delpc(1,    1) = delpc(1,    1) + vort(1,    0)
          if (se_corner) delpc(npx,  1) = delpc(npx,  1) + vort(npx,  0)
          if (ne_corner) delpc(npx,npy) = delpc(npx,npy) - vort(npx,npy)
          if (nw_corner) delpc(1,  npy) = delpc(1,  npy) - vort(1,  npy)

          damp = dddmp*da_min_c

          do j=js,je+1
             do i=is,ie+1
!               ke(i,j) = ke(i,j) - dddmp*delpc(i,j)  ! divg defined at corners
                ke(i,j) = ke(i,j) - damp*rarea_c(i,j)*delpc(i,j)
             enddo
          enddo


! Calc Vorticity
! Convert winds to circulation elements:
       do j=jsd,jed+1
          do i=isd,ied
             u(i,j) = u(i,j)*dx(i,j)
          enddo
       enddo
       do j=jsd,jed
          do i=isd,ied+1
             v(i,j) = v(i,j)*dy(i,j)
          enddo
       enddo

       do j=jsd,jed
          do i=isd,ied
#ifdef NON_ROTATING
             vort(i,j) = rarea(i,j)*(u(i,j)-u(i,j+1)-v(i,j)+v(i+1,j))
#else
             vort(i,j) = f0(i,j) + rarea(i,j)*(u(i,j)-u(i,j+1)-v(i,j)+v(i+1,j))
#endif
          enddo
       enddo

       call fv_tp_2d(vort, crx_adv, cry_adv, npx-1, npy-1, hord_vt, &
                     fx, fy, xfx_adv,yfx_adv, area, uniform_ppm)

       do j=js,je+1
          do i=is,ie
             u(i,j) = u(i,j) + ke(i,j) - ke(i+1,j) + fy(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j) = v(i,j) + ke(i,j) - ke(i,j+1) - fx(i,j)
          enddo
       enddo

#ifdef SW_DYNAMICS
      endif ! test_case
#endif

   end subroutine d_sw

 subroutine xtp_u(c, u, v, flux, xfx, iord)

 real, INTENT(IN)  ::   u(isd:ied,jsd:jed+1)
 real, INTENT(IN)  ::   v(isd:ied+1,jsd:jed) ! u-wind
 real, INTENT(IN)  ::   c(is:ie+1,js:je+1)   !  Courant   N (like FLUX)
 real, INTENT(in)  :: xfx(is:ie+1,js:je+1)
 real, INTENT(out):: flux(is:ie+1,js:je+1)
 integer, INTENT(IN) ::  iord
! Local
 real al(is-1:ie+2), dm(is-2:ie+2)
 real bl(is-1:ie+1)
 real br(is-1:ie+1)
 real dq(is-3:ie+2)
 real dl, dr, xt, pmp, lac, dqt, q0
 integer i, j

 if (iord<=4) then
     do j=js,je+1

        do i=is-2,ie+2
           xt = 0.25*(u(i+1,j) - u(i-1,j))
           dm(i) = sign(min(abs(xt), max(u(i-1,j), u(i,j), u(i+1,j)) - u(i,j),  &
                            u(i,j) - min(u(i-1,j), u(i,j), u(i+1,j))), xt)
        enddo

        if ( j==1 .or. j==npy ) then      ! top & bottom edges
            if(  is   ==1 )    dm(1 ) = 0.
            if( (ie+1)==npx )  dm(ie) = 0.
        endif

        do i=is-1,ie+2
           al(i) = 0.5*(u(i-1,j)+u(i,j)) + r3*(dm(i-1) - dm(i))
        enddo

        do i=is,ie+1
          if( c(i,j)>0. ) then
             xt = 2.*dm(i-1)
             dl = sign(min(abs(xt), abs(al(i-1)-u(i-1,j))), xt)
             dr = sign(min(abs(xt), abs(al(i  )-u(i-1,j))), xt)
             flux(i,j) = u(i-1,j) + (1.-c(i,j))*(dr + c(i,j)*(dl-dr))
          else
             xt = 2.*dm(i)
             dl = sign(min(abs(xt), abs(al(i  )-u(i,j))), xt)
             dr = sign(min(abs(xt), abs(al(i+1)-u(i,j))), xt)
             flux(i,j) = u(i,j) - (1.+c(i,j))*(dl + c(i,j)*(dl-dr))
          endif 
        enddo
     enddo
 elseif ( iord==5 .or. iord==6 ) then
     do j=js,je+1
        do i=is-3,ie+2
           dq(i) = u(i+1,j) - u(i,j)
        enddo
        do i=is-2,ie+2
              xt = 0.25*(dq(i-1) + dq(i))
           dm(i) = sign(min(abs(xt), abs(dq(i-1)), abs(dq(i))), xt)
        enddo
        do i=max(3,is-1),min(npx-2,ie+2)
           al(i) = 0.5*(u(i-1,j)+u(i,j)) + r3*(dm(i-1) - dm(i))
        enddo

! Perturbation form:
        do i=max(3,is-1),min(npx-3,ie+1)
           pmp = -2.*dq(i)
           lac = pmp + 1.5*dq(i+1)
           bl(i) = min(max(0., pmp, lac), max(al(i  )-u(i,j), min(0.,pmp, lac)))
           pmp = 2.*dq(i-1)
           lac = pmp - 1.5*dq(i-2)
           br(i) = min(max(0., pmp, lac), max(al(i+1)-u(i,j), min(0.,pmp, lac)))
        enddo

!--------------
! fix the edges
!--------------
        if ( is==1 ) then
             br(2) = al(3) - u(2,j)
                xt = 3./14.*u(1,j) + 11./14.*u( 2,j) - 4./7.*dm( 2)
             bl(2) = xt - u(2,j)

           if( j==1 .or. j==npy ) then
              bl(0) = 0.
              br(0) = 0.
              bl(1) = 0.
              br(1) = 0.
              if(iord==5) call std_ppm(1, u(2,j), bl(2), br(2))
           else
              br(1) = xt - u(1,j)
              bl(0) = 11./14.*(u(-1,j)-u(0,j)) + 4./7.*dm(-1)
              xt = 27./28.*(u(0,j)+u(1,j)) - 13./28.*(u(-1,j)+u(2,j))   &
                  + 3./7.*(dm(2)-dm(-1))
#ifdef TEST_1E
              dl = 0.5*(v(1,j-1)+v(1,j)) * cosa(1,j) 
              bl(1) = xt + dl - u(1,j)
              br(0) = xt - dl - u(0,j)
#else
              bl(1) = xt - u(1,j)
              br(0) = xt - u(0,j)
#endif
              if(iord==5) call std_ppm(3, u(0,j), bl(0), br(0))
           endif
        endif
        if ( (ie+1)==npx ) then
             bl(npx-2) = al(npx-2) - u(npx-2,j)
                    xt = 3./14.*u(npx-1,j) + 11./14.*u(npx-2,j) + 4./7.*dm(npx-2)
             br(npx-2) = xt - u(npx-2,j)
           if( j==1 .or. j==npy ) then
             bl(npx-1) = 0.
             br(npx-1) = 0.
             bl(npx  ) = 0.
             br(npx  ) = 0.
             if(iord==5) call std_ppm(1, u(npx-2,j), bl(npx-2), br(npx-2))
           else
             bl(npx-1) = xt - u(npx-1,j)
             br(npx) = 11./14.*(u(npx+1,j)-u(npx,j)) - 4./7.*dm(npx+1)
             xt = 27./28.*(u(npx-1,j)+u(npx,j)) - 13./28.*(u(npx-2,j)+u(npx+1,j))   &
                 + 3./7.*(dm(npx+1)-dm(npx-2))
#ifdef TEST_1E
             dl = 0.5*(v(npx,j-1)+v(npx,j)) * cosa(npx,j) 
             br(npx-1) = xt + dl - u(npx-1,j)
             bl(npx  ) = xt - dl - u(npx  ,j)
#else
             br(npx-1) = xt - u(npx-1,j)
             bl(npx  ) = xt - u(npx  ,j)
#endif
             if(iord==5) call std_ppm(3, u(npx-2,j), bl(npx-2), br(npx-2))
           endif
        endif

        do i=is,ie+1
           if( c(i,j)>0. ) then
              flux(i,j) = u(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
           else
              flux(i,j) = u(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
           endif
        enddo
     enddo
 else
     do j=js,je+1
        do i=is-2,ie+2
              xt = 0.25*(u(i+1,j) - u(i-1,j))
           dm(i) = sign(min(abs(xt), max(u(i-1,j), u(i,j), u(i+1,j)) - u(i,j),  &
                            u(i,j) - min(u(i-1,j), u(i,j), u(i+1,j))), xt)
        enddo

        do i=max(3,is-1),min(npx-2,ie+2)
           al(i) = 0.5*(u(i-1,j)+u(i,j)) + r3*(dm(i-1) - dm(i))
        enddo

        do i=max(3,is-1),min(npx-3,ie+1)
              xt = 2.*dm(i)
           bl(i) = sign(min(abs(xt), abs(al(i)  -u(i,j))), xt)
           br(i) = sign(min(abs(xt), abs(al(i+1)-u(i,j))), xt)
        enddo

!--------------
! fix the edges
!--------------
             if ( is==1 ) then
! Inside:
                   dqt = u(1,j) - u(2,j)
                   bl(1) =  6./7.*dm(2) + 13./14.*dqt
                   br(1) = -4./7.*dm(2) - 11./14.*dqt
                      xt = 0.5*(br(1) - bl(1))
                  dm(1) = sign(min(abs(xt), abs(dqt)), xt)
                     xt = 2.*dm(1)
                  bl(1) = sign( min(abs(xt), abs(bl(1))), xt )
                  br(1) = sign( min(abs(xt), abs(br(1))), xt )
                  bl(2) = 0.5*dqt + r3*(dm(1) - dm(2))
                     xt = 2.*dm(2)
                  bl(2) = sign(min(abs(xt), abs(bl(2))), xt)
                  br(2) = sign(min(abs(xt), abs(al(3)-u(2,j))), xt)
! Outside: (mirrors what's done inside)
                   dqt = u(0,j) - u(-1,j)
                   br(0) = -6./7.*dm(-1) + 13./14.*dqt
                   bl(0) =  4./7.*dm(-1) - 11./14.*dqt
                      xt = 0.5*(br(0) - bl(0))
                  dm(0) = sign(min(abs(xt), abs(dqt)), xt)
                     xt = 2.*dm(0)
                  bl(0) = sign( min(abs(xt), abs(bl(0))), xt )
                  br(0) = sign( min(abs(xt), abs(br(0))), xt )
! Average (u+,u-) at edges:
!                 xt  = 0.5*(bl(1)+br(0)+u(0,j)+u(1,j))
!                 dqt = 0.5*(v(1,j-1)+v(1,j))*cosa(1,j) 
!                 bl(1) =  (xt + dqt - u(1,j))
!                 br(0) =   xt - dqt - u(0,j)
             endif
             if ( (ie+1)==npx ) then
! Inside:
                   dqt = u(npx-1,j) - u(npx-2,j)
                   br(npx-1) = -6./7.*dm(npx-2) + 13./14.*dqt
                   bl(npx-1) =  4./7.*dm(npx-2) - 11./14.*dqt
                      xt = 0.5*(br(npx-1) - bl(npx-1))
                  dm(npx-1) = sign(min(abs(xt), abs(dqt)), xt)
                     xt = 2.*dm(npx-1)
                  bl(npx-1) = sign( min(abs(xt), abs(bl(npx-1))), xt )
                  br(npx-1) = sign( min(abs(xt), abs(br(npx-1))), xt )
                  br(npx-2) = 0.5*dqt + r3*(dm(npx-2) - dm(npx-1))
                         xt = 2.*dm(npx-2)
                  br(npx-2) = sign(min(abs(xt), abs(br(npx-2))), xt)
                  bl(npx-2) = sign(min(abs(xt), abs(al(npx-2)-u(npx-2,j))), xt)
! Outside: (mirrors what's done inside)
                   dqt = u(npx,j) - u(npx+1,j)
                   bl(npx) =  6./7.*dm(npx+1) + 13./14.*dqt
                   br(npx) = -4./7.*dm(npx+1) - 11./14.*dqt
                        xt = 0.5*(br(npx) - bl(npx))
                  dm(npx) = sign(min(abs(xt), abs(dqt)), xt)
                     xt = 2.*dm(npx)
                  bl(npx) = sign( min(abs(xt), abs(bl(npx))), xt )
                  br(npx) = sign( min(abs(xt), abs(br(npx))), xt )
! Average (u+,u-) at Eastern edge:
!                 xt  = 0.5*(br(npx-1)+bl(npx)+u(npx-1,j)+u(npx,j))
!                 dqt = 0.5*(v(npx,j-1)+v(npx,j))*cosa(npx,j) 
!                 br(npx-1) =   xt + dqt - u(npx-1,j)
!                 bl(npx  ) =  (xt - dqt - u(npx,  j))
        endif

        do i=is,ie+1
           if( c(i,j)>0. ) then
              flux(i,j) = u(i-1,j) + (1.-c(i,j))*(br(i-1)+c(i,j)*(bl(i-1)-br(i-1)))
           else
              flux(i,j) = u(i,  j) - (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )-br(i  )))
           endif
        enddo
     enddo
 endif

   do j=js,je+1
      do i=is,ie+1
         flux(i,j) = flux(i,j)*xfx(i,j)
      enddo
   enddo

 end subroutine xtp_u


 subroutine ytp_v(c, u, v, flux, yfx, jord)
 integer, intent(IN):: jord
 real, INTENT(IN)  ::   u(isd:ied,jsd:jed+1)
 real, INTENT(IN)  ::   v(isd:ied+1,jsd:jed) ! u-wind
 real, INTENT(IN) ::    c(is:ie+1,js:je+1)   !  Courant   N (like FLUX)
 real, INTENT(IN) ::  yfx(is:ie+1,js:je+1)
 real, INTENT(OUT):: flux(is:ie+1,js:je+1)
! Local:
 real dm(is:ie+1,js-2:je+2)
 real al(is:ie+1,js-1:je+2)
 real bl(is:ie+1,js-1:je+1)
 real br(is:ie+1,js-1:je+1)
 real dq(is:ie+1,js-3:je+2)
 real xt, dl, dr, pmp, lac, dqt, q0
 integer i, j

 if (jord<=4) then

     do j=js-2,je+2
        do i=is,ie+1
                xt = 0.25*(v(i,j+1) - v(i,j-1))
           dm(i,j) = sign(min(abs(xt), max(v(i,j-1), v(i,j), v(i,j+1)) - v(i,j),   &
                              v(i,j) - min(v(i,j-1), v(i,j), v(i,j+1))), xt)
        enddo
     enddo

     if ( is==1 ) then
         if (  js   ==1 )   dm(1, 1) = 0.
         if ( (je+1)==npy ) dm(1,je) = 0.
     endif

     if ( (ie+1)==npx ) then
         if (  js   ==1   )  dm(npx, 1) = 0.
         if ( (je+1)==npy )  dm(npx,je) = 0.
     endif


   do j=js-1,je+2
      do i=is,ie+1
         al(i,j) = 0.5*(v(i,j-1)+v(i,j)) + r3*(dm(i,j-1) - dm(i,j))
      enddo
   enddo

   do j=js,je+1
      do i=is,ie+1
         if(c(i,j)>0.) then
            xt = 2.*dm(i,j-1)
            dl = sign(min(abs(xt), abs(al(i,j-1)-v(i,j-1))), xt)
            dr = sign(min(abs(xt), abs(al(i,j)-v(i,j-1))),   xt)
            flux(i,j) = v(i,j-1) + (1.-c(i,j))*(dr + c(i,j)*(dl-dr))
         else
            xt = 2.*dm(i,j)
            dl = sign(min(abs(xt), abs(al(i,j)-v(i,j))),   xt)
            dr = sign(min(abs(xt), abs(al(i,j+1)-v(i,j))), xt)
            flux(i,j) = v(i,j) - (1.+c(i,j))*(dl + c(i,j)*(dl-dr))
         endif
      enddo
   enddo
 elseif (jord==5 .or. jord==6) then
! PPM with Hunyh's 2nd constraint
   do j=js-3, je+2
      do i=is,ie+1
         dq(i,j) = v(i,j+1) - v(i,j)
      enddo
   enddo

   do j=js-2,je+2
      do i=is,ie+1
         xt = 0.25*(dq(i,j-1) + dq(i,j))
         dm(i,j) = sign(min(abs(xt), abs(dq(i,j-1)), abs(dq(i,j))), xt)
      enddo
   enddo

   do j=max(3,js-1),min(npy-2,je+2)
      do i=is,ie+1
         al(i,j) = 0.5*(v(i,j-1)+v(i,j)) + r3*(dm(i,j-1) - dm(i,j))
      enddo
   enddo

   do j=max(3,js-1),min(npy-3,je+1)
      do i=is,ie+1
            pmp = -2.*dq(i,j) 
            lac = pmp + 1.5*dq(i,j+1)
            bl(i,j) = min(max(0.,pmp,lac), max(al(i,j)-v(i,j), min(0.,pmp,lac)))
            pmp = 2.*dq(i,j-1)
            lac = pmp - 1.5*dq(i,j-2)
            br(i,j) = min(max(0.,pmp,lac), max(al(i,j+1)-v(i,j), min(0.,pmp,lac)))
      enddo
   enddo

   if( js==1 ) then
!      do i=max(2,is),min(npx-1,ie+1)
       do i=is,ie+1
          br(i,2) = al(i,3) - v(i,2)
               xt = 3./14.*v(i,1) + 11./14.*v(i,2) - 4./7.*dm(i,2)
          bl(i,2) = xt - v(i,2)
          br(i,1) = xt - v(i,1)
          bl(i,0) = 11./14.*(v(i,-1)-v(i,0)) + 4./7.*dm(i,-1)
             xt = 27./28.*(v(i,0)+v(i,1)) - 13./28.*(v(i,-1)+v(i,2))   &
                 + 3./7.*(dm(i,2)-dm(i,-1))
#ifdef TEST_1E
          dl = 0.5*(u(i-1,1)+u(i,1)) * cosa(i,1) 
          bl(i,1) = xt + dl - v(i,1)
          br(i,0) = xt - dl - v(i,0)
#else
          bl(i,1) = xt - v(i,1)
          br(i,0) = xt - v(i,0)
#endif
       enddo
       if ( jord==5 ) then
          do j=0,2
             call std_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j))
          enddo
       endif

       if ( is==1 ) then
          do j=0,1
             bl(1,j) = 0.
             br(1,j) = 0.
          enddo
       endif
       if ( (ie+1)==npx ) then
          do j=0,1
             bl(npx,j) = 0.
             br(npx,j) = 0.
          enddo
       endif
   endif 

   if( (je+1)==npy ) then
!      do i=max(2,is),min(npx-1,ie+1)
       do i=is,ie+1
          bl(i,npy-2) = al(i,npy-2) - v(i,npy-2)
                   xt = 3./14.*v(i,npy-1) + 11./14.*v(i,npy-2) + 4./7.*dm(i,npy-2)
          br(i,npy-2) = xt - v(i,npy-2)

          bl(i,npy-1) = xt - v(i,npy-1)
          br(i,npy) = 11./14.*(v(i,npy+1)-v(i,npy)) - 4./7.*dm(i,npy+1)
          xt = 27./28.*(v(i,npy-1)+v(i,npy)) - 13./28.*(v(i,npy-2)+v(i,npy+1))   &
              + 3./7.*(dm(i,npy+1)-dm(i,npy-2))
#ifdef TEST_1E
          dl = 0.5*(u(i-1,npy)+u(i,npy)) * cosa(i,npy) 
          br(i,npy-1) = xt + dl - v(i,npy-1)
          bl(i,npy  ) = xt - dl - v(i,npy)
#else
          br(i,npy-1) = xt - v(i,npy-1)
          bl(i,npy  ) = xt - v(i,npy)
#endif
       enddo
       if ( jord==5 ) then
          do j=npy-2,npy
             call std_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j))
          enddo
       endif

        if ( is==1 ) then
          do j=npy-1,npy
             bl(1,j) = 0.
             br(1,j) = 0.
          enddo
        endif
        if ( (ie+1)==npx ) then
          do j=npy-1,npy
             bl(npx,j) = 0.
             br(npx,j) = 0.
          enddo
        endif
   endif

   do j=js,je+1
      do i=is,ie+1
         if(c(i,j)>0.) then
            flux(i,j) = v(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = v(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo

 else
   do j=js-2,je+2
      do i=is,ie+1
              xt = 0.25*(v(i,j+1) - v(i,j-1))
         dm(i,j) = sign(min(abs(xt), max(v(i,j-1), v(i,j), v(i,j+1)) - v(i,j),   &
                            v(i,j) - min(v(i,j-1), v(i,j), v(i,j+1))), xt)
      enddo
   enddo

   do j=max(3,js-1),min(npy-2,je+2)
      do i=is,ie+1
         al(i,j) = 0.5*(v(i,j-1)+v(i,j)) + r3*(dm(i,j-1) - dm(i,j))
      enddo
   enddo

   do j=max(3,js-1),min(npy-3,je+1)
      do i=is,ie+1
              xt = 2.*dm(i,j)
         bl(i,j) = sign(min(abs(xt), abs(al(i,j  )-v(i,j))), xt)
         br(i,j) = sign(min(abs(xt), abs(al(i,j+1)-v(i,j))), xt)
      enddo
   enddo

      if( js==1 ) then
        do i=is,ie+1
! Inside:
          dqt = v(i,1) - v(i,2)
          bl(i,1) =  6./7.*dm(i,2) + 13./14.*dqt
          br(i,1) = -4./7.*dm(i,2) - 11./14.*dqt
               xt = 0.5*(br(i,1) - bl(i,1))
          dm(i,1) = sign(min(abs(xt), abs(dqt)), xt)
!----------------------------------------------------
!         br(i,1) = -0.5*dqt + r3*(dm(i,1) - dm(i,2))
!----------------------------------------------------
               xt = 2.*dm(i,1)
          bl(i,1) = sign( min(abs(xt), abs(bl(i,1))), xt )
          br(i,1) = sign( min(abs(xt), abs(br(i,1))), xt )
          bl(i,2) = 0.5*dqt + r3*(dm(i,1) - dm(i,2))
               xt = 2.*dm(i,2)
          bl(i,2) = sign(min(abs(xt), abs(bl(i,2))), xt)
          br(i,2) = sign(min(abs(xt), abs(al(i,3)-v(i,2))), xt)
! Outside: (mirrors what's done inside)
              dqt = v(i,0) - v(i,-1)
          br(i,0) = -6./7.*dm(i,-1) + 13./14.*dqt
          bl(i,0) =  4./7.*dm(i,-1) - 11./14.*dqt
               xt = 0.5*(br(i,0) - bl(i,0))
          dm(i,0) = sign(min(abs(xt), abs(dqt)), xt)
               xt = 2.*dm(i,0)
          bl(i,0) = sign( min(abs(xt), abs(bl(i,0))), xt )
          br(i,0) = sign( min(abs(xt), abs(br(i,0))), xt )
!-------------------------
! Average (v+,v-) at south edge (j=1):
!             xt  = 0.5*(bl(i,1)+br(i,0)+v(i,0)+v(i,1))
!             dqt = 0.5*(u(i-1,1)+u(i,1))*cosa(i,1) 
!         bl(i,1) =  (xt + dqt - v(i,1))
!         br(i,0) =   xt - dqt - v(i,0)
        enddo
      endif 
      if( (je+1)==npy ) then
        do i=is,ie+1
! Inside:
                   dqt = v(i,npy-1) - v(i,npy-2)
           br(i,npy-1) = -6./7.*dm(i,npy-2) + 13./14.*dqt
           bl(i,npy-1) =  4./7.*dm(i,npy-2) - 11./14.*dqt
                    xt = 0.5*(br(i,npy-1) - bl(i,npy-1))
           dm(i,npy-1) = sign(min(abs(xt), abs(dqt)), xt)
                    xt = 2.*dm(i,npy-1)
           bl(i,npy-1) = sign( min(abs(xt), abs(bl(i,npy-1))), xt )
           br(i,npy-1) = sign( min(abs(xt), abs(br(i,npy-1))), xt )
           br(i,npy-2) = 0.5*dqt + r3*(dm(i,npy-2) - dm(i,npy-1))
                    xt = 2.*dm(i,npy-2)
           br(i,npy-2) = sign(min(abs(xt), abs(br(i,npy-2))), xt)
           bl(i,npy-2) = sign(min(abs(xt), abs(al(i,npy-2)-v(i,npy-2))), xt)
! Outside: (mirrors what's done inside)
                 dqt = v(i,npy) - v(i,npy+1)
           bl(i,npy) =  6./7.*dm(i,npy+1) + 13./14.*dqt
           br(i,npy) = -4./7.*dm(i,npy+1) - 11./14.*dqt
                  xt = 0.5*(br(i,npy) - bl(i,npy))
           dm(i,npy) = sign(min(abs(xt), abs(dqt)), xt)
                  xt = 2.*dm(i,npy)
           bl(i,npy) = sign( min(abs(xt), abs(bl(i,npy))), xt )
           br(i,npy) = sign( min(abs(xt), abs(br(i,npy))), xt )
!-------------------------
! Average (v+,v-) at north edge (j=npy):
!             xt  = 0.5*(br(i,npy-1)+bl(i,npy)+v(i,npy-1)+v(i,npy))
!             dqt = 0.5*(u(i-1,npy)+u(i,npy))*cosa(i,npy) 
!         br(i,npy-1) =   xt + dqt - v(i,npy-1)
!         bl(i,npy  ) =  (xt - dqt - v(i,npy))
        enddo
      endif

   do j=js,je+1
      do i=is,ie+1
         if(c(i,j)>0.) then
            flux(i,j) = v(i,j-1) + (1.-c(i,j))*(br(i,j-1)+c(i,j)*(bl(i,j-1)-br(i,j-1)))
         else
            flux(i,j) = v(i,j  ) - (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )-br(i,j  )))
         endif
      enddo
   enddo
 endif

   do j=js,je+1
      do i=is,ie+1
         flux(i,j) = flux(i,j)*yfx(i,j)
      enddo
   enddo

 end subroutine ytp_v



 subroutine d2a2c_vect(u, v, ua, va, uc, vc, ut, vt)
  real, intent(in) ::  u(isd:ied,jsd:jed+1)
  real, intent(in) ::  v(isd:ied+1,jsd:jed)
  real, intent(out), dimension(isd:ied+1,jsd:jed  ):: uc
  real, intent(out), dimension(isd:ied  ,jsd:jed+1):: vc
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: ua, va, ut, vt
! Local 
!----------------------------------------------
! 4-pt Lagrange interpolation
#ifdef LAGR_INTP
  real, parameter:: a1 =  0.5625
  real, parameter:: a2 = -0.0625
#else
! PPM volume mean form:
  real, parameter:: a1 =   7./12.     ! 0.58333333
  real, parameter:: a2 =  -1./12.
#endif
!----------------------------------------------
! 3-pt formular:
  real, parameter:: b1 = -0.125
  real, parameter:: b2 =  0.75
  real, parameter:: b3 =  0.375
!----------------------------------------------
  real, dimension(isd:ied,jsd:jed):: utmp, vtmp
  integer i, j, ifirst, ilast

  do j=jsd,jed
     do i=isd,ied+1
        uc(i,j) = v(i,j)*dy(i,j)
     enddo
  enddo

  do j=jsd,jed+1
     do i=isd,ied
        vc(i,j) = u(i,j)*dx(i,j)
     enddo
  enddo

! D --> A
! Co-variant to Co-variant "vorticity-conserving" remapping
  do j=jsd,jed
     do i=isd,ied
        utmp(i,j) = 0.5*(vc(i,j) + vc(i,j+1)) * rdxa(i,j)
        vtmp(i,j) = 0.5*(uc(i,j) + uc(i+1,j)) * rdya(i,j)
        ua(i,j) = (utmp(i,j)-vtmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
        va(i,j) = (vtmp(i,j)-utmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
     enddo
  enddo

! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
     if( sw_corner ) then
        utmp( 0, 0) = -vtmp(0, 1)
        utmp(-1, 0) = -vtmp(0, 2)
!       utmp(-2, 0) = -vtmp(0, 3)
     endif
     if( se_corner ) then
        utmp(npx,  0) = vtmp(npx, 1)
        utmp(npx+1,0) = vtmp(npx, 2)
!       utmp(npx+2,0) = vtmp(npx, 3)
     endif
     if( ne_corner ) then
         utmp(npx,  npy) = -vtmp(npx,je  )
         utmp(npx+1,npy) = -vtmp(npx,je-1)
!        utmp(npx+2,npy) = -vtmp(npx,je-2)
     endif
     if( nw_corner ) then
         utmp( 0,npy) = vtmp(0,je  )
         utmp(-1,npy) = vtmp(0,je-1)
!        utmp(-2,npy) = vtmp(0,je-2)
     endif

     ifirst = max(3,    is-1)
     ilast  = min(npx-2,ie+2)

! 4th order interpolation for interior points:
     do j=js-1,je+1
        do i=ifirst,ilast
           uc(i,j) = a1*(utmp(i-1,j) + utmp(i,  j)) +  &
                     a2*(utmp(i-2,j) + utmp(i+1,j))
           ut(i,j) = (uc(i,j) - v(i,j)*cosa_u(i,j))*rsin_u(i,j)
        enddo
     enddo

     if( is==1 ) then
        do j=js-1,je+1
! Special case: i=1
! 2-pt extrapolation --------------------------------------------------------
           uc(1,j) = 0.25*( -utmp(-1,j) + 3.*(utmp(0,j)+utmp(1,j))     &
                            -utmp( 2,j) ) * rsin_u(1,j)
! 3-pt extrapolation --------------------------------------------------------
!          uc(1,j) = ( 24.*(utmp(0,j)+utmp(1,j)) - 13.*(utmp(-1,j)+utmp(2,j))    &
!                     + 3.*(utmp(-2,j)+utmp(3,j)) ) * rsin_u(1,j) / 28.
           ut(1,j) = uc(1,j) * rsin_u(1,j)
! 3-pt formular:
           uc(0,j) = b1*utmp(-2,j) + b2*utmp(-1,j) + b3*utmp(0,j) 
           ut(0,j) = (uc(0,j) - v(0,j)*cosa_u(0,j))*rsin_u(0,j)
           uc(2,j) = b3*utmp(1,j) + b2*utmp(2,j) + b1*utmp(3,j)
           ut(2,j) = (uc(2,j) - v(2,j)*cosa_u(2,j))*rsin_u(2,j)
        enddo
     endif
     if( (ie+1)==npx ) then
        do j=js-1,je+1
! 2-pt extrapolation --------------------------------------------------------
           uc(npx,j) = 0.25*(-utmp(npx-2,j) + 3.*(utmp(ie,j)+utmp(npx,j))  &
                             -utmp(npx+1,j) ) * rsin_u(npx,j)
! 3-pt extrapolation --------------------------------------------------------
!          uc(npx,j) = ( 24.*(utmp(npx-1,j)+utmp(npx,j))-13.*(utmp(npx-2,j)+utmp(npx+1,j))    &
!                       + 3.*(utmp(npx-3,j)+utmp(npx+2,j)) ) * rsin_u(npx,j) / 28.
           ut(npx,j) = uc(npx,j) * rsin_u(npx,j)
! 3-pt formular:
           uc(npx-1,j) = b1*utmp(npx-3,j) + b2*utmp(npx-2,j) + b3*utmp(npx-1,j) 
           ut(npx-1,j) = (uc(npx-1,j) - v(npx-1,j)*cosa_u(npx-1,j))*rsin_u(npx-1,j)
           uc(npx+1,j) = b3*utmp(npx,j) + b2*utmp(npx+1,j) + b1*utmp(npx+2,j) 
           ut(npx+1,j) = (uc(npx+1,j) - v(npx+1,j)*cosa_u(npx+1,j))*rsin_u(npx+1,j)
        enddo
     endif

!------
! Ydir:
!------
     if( sw_corner ) then
         vtmp(0, 0) = -utmp(1,0)
         vtmp(0,-1) = -utmp(2,0)
!        vtmp(0,-2) = -utmp(3,0)
     endif
     if( nw_corner ) then
         vtmp(0, npy  ) =  utmp(1,npy)
         vtmp(0, npy+1) =  utmp(2,npy)
!        vtmp(0, npy+2) =  utmp(3,npy)
     endif
     if( se_corner ) then
         vtmp(npx, 0) =  utmp(ie,  0)
         vtmp(npx,-1) =  utmp(ie-1,0)
!        vtmp(npx,-2) =  utmp(ie-2,0)
     endif
     if( ne_corner ) then
         vtmp(npx,npy  ) = -utmp(ie,  npy)
         vtmp(npx,npy+1) = -utmp(ie-1,npy)
!        vtmp(npx,npy+2) = -utmp(ie-2,npy)
     endif

     do j=js-1,je+2
      if ( j==1 ) then
        do i=is-1,ie+1
! 2-pt extrapolation --------------------------------------------------------
           vc(i,1) = 0.25 * ( -vtmp(i,-1) +  3.*(vtmp(i,0)+vtmp(i,1))       &
                              -vtmp(i,2) ) * rsin_v(i,1)
! 3-pt extrapolation --------------------------------------------------------
!          vc(i,1) = (24.*(vtmp(i,0)+vtmp(i,1))-13.*(vtmp(i,-1)+vtmp(i,2))  &
!                    + 3.*(vtmp(i,-2)+vtmp(i,3)) ) * rsin_v(i,1) / 28.
           vt(i,1) = vc(i,1) * rsin_v(i,1)
        enddo
      elseif ( j==0 .or.  j==(npy-1) ) then
        do i=is-1,ie+1
           vc(i,j) = b1*vtmp(i,j-2) + b2*vtmp(i,j-1) + b3*vtmp(i,j)
           vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
        enddo
      elseif ( j==2 .or. j==(npy+1) ) then
        do i=is-1,ie+1
           vc(i,j) = b3*vtmp(i,j-1) + b2*vtmp(i,j) + b1*vtmp(i,j+1)
           vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
        enddo
      elseif ( j==npy ) then
        do i=is-1,ie+1
! 2-pt extrapolation --------------------------------------------------------
           vc(i,npy) = 0.25*( -vtmp(i,npy-2)+ 3.*(vtmp(i,je)+vtmp(i,npy))   &
                              -vtmp(i,npy+1) )*rsin_v(i,npy)
! 3-pt extrapolation --------------------------------------------------------
!          vc(i,npy) = (24.*(vtmp(i,npy-1)+vtmp(i,npy))-13.*(vtmp(i,npy-2)+vtmp(i,npy+1))  &
!                      + 3.*(vtmp(i,npy-3)+vtmp(i,npy+2)))*rsin_v(i,npy) / 28.
           vt(i,npy) = vc(i,npy) * rsin_v(i,npy)
        enddo
      else
! 4th order interpolation for interior points:
        do i=is-1,ie+1
           vc(i,j) = a1*(vtmp(i,j-1)+vtmp(i,j  )) +   &
                     a2*(vtmp(i,j-2)+vtmp(i,j+1))
           vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
        enddo
      endif
     enddo

 end subroutine d2a2c_vect
 


 subroutine d2a2c_vect_v0( u, v, ua, va, uc, vc, ut, vt )
  real, intent(in) ::  u(isd:ied,jsd:jed+1)
  real, intent(in) ::  v(isd:ied+1,jsd:jed)
  real, intent(out), dimension(isd:ied+1,jsd:jed  ):: uc
  real, intent(out), dimension(isd:ied  ,jsd:jed+1):: vc
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: ua, va, ut, vt
! Local 
    real, dimension(is-2:ie+2,js-2:je+2):: wk
#ifdef VECT_DEV
    real, dimension(isd:ied,jsd:jed):: utmp, vtmp
#else
    real :: utmp, vtmp
#endif
    real, dimension(isd:ied,jsd:jed):: v1, v2, v3
    real v1d(3,-3:npx+4), w3(3,-3:npx+4)
    real vw1, vw2, vw3
    real vs1, vs2, vs3
    integer k, im, jm, im2, jm2
    integer i, j


! needs only ut[is-1:ie+2,js-1:je+1], vt[is-1:ie+1,js-1:je+2]
     do j=js-2,je+2
        do i=is-2,ie+3
           ut(i,j) = v(i,j)*dy(i,j)
        enddo
     enddo

     do j=js-2,je+3
        do i=is-2,ie+2
           vt(i,j) = u(i,j)*dx(i,j)
        enddo
     enddo

! D --> A
     do j=js-2,je+2
        do i=is-2,ie+2
! Co-variant to Co-variant "vorticity-conserving" interpolation
#ifdef VECT_DEV
           utmp(i,j) = 0.5*(vt(i,j) + vt(i,j+1)) * rdxa(i,j)
           vtmp(i,j) = 0.5*(ut(i,j) + ut(i+1,j)) * rdya(i,j)
           ua(i,j) = (utmp(i,j)-vtmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
           va(i,j) = (vtmp(i,j)-utmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
#else
           utmp = 0.5*(vt(i,j) + vt(i,j+1)) * rdxa(i,j)
           vtmp = 0.5*(ut(i,j) + ut(i+1,j)) * rdya(i,j)
           ua(i,j) = (utmp-vtmp*cosa_s(i,j)) * rsin2(i,j)
           va(i,j) = (vtmp-utmp*cosa_s(i,j)) * rsin2(i,j)
#endif
        enddo
     enddo

! A -> C
!--------------------------------------------
! Divergence conserving interp to cell walls
!--------------------------------------------
     do j=js-1,je+1
        do i=is-2,ie+2
           wk(i,j) = ua(i,j)*dya(i,j)*sina_s(i,j)
        enddo
     enddo
     do j=js-1,je+1
        do i=is-1,ie+2
           ut(i,j) = 0.5*(wk(i-1,j)+wk(i,j)) / (dy(i,j)*sina_u(i,j))
           uc(i,j) = ut(i,j) + 0.5*(va(i-1,j)*cosa_s(i-1,j)+va(i,j)*cosa_s(i,j))
        enddo
     enddo

     do j=js-2,je+2
        do i=is-1,ie+1
           wk(i,j) = va(i,j)*dxa(i,j)*sina_s(i,j)
        enddo
     enddo
     do j=js-1,je+2
        do i=is-1,ie+1
           vt(i,j) = 0.5*(wk(i,j-1)+wk(i,j)) / (dx(i,j)*sina_v(i,j))
           vc(i,j) = vt(i,j) + 0.5*(ua(i,j-1)*cosa_s(i,j-1)+ua(i,j)*cosa_s(i,j))
        enddo
     enddo

!--------------
! Fix the edges
!--------------
! V = ua * e1 + va * e2
     do j=js-2,je+2
     if ( is==1 ) then
        do i=0,1
        v1(i,j) = ua(i,j)*ec1(1,i,j) + va(i,j)*ec2(1,i,j)
        v2(i,j) = ua(i,j)*ec1(2,i,j) + va(i,j)*ec2(2,i,j)
        v3(i,j) = ua(i,j)*ec1(3,i,j) + va(i,j)*ec2(3,i,j)
        enddo
     endif
     if ( (ie+1)==npx ) then
        do i=npx-1,npx
        v1(i,j) = ua(i,j)*ec1(1,i,j) + va(i,j)*ec2(1,i,j)
        v2(i,j) = ua(i,j)*ec1(2,i,j) + va(i,j)*ec2(2,i,j)
        v3(i,j) = ua(i,j)*ec1(3,i,j) + va(i,j)*ec2(3,i,j)
        enddo
     endif
     enddo

     if ( js==1 ) then
     do j=0,1
        do i=is-2,ie+2
           v1(i,j) = ua(i,j)*ec1(1,i,j) + va(i,j)*ec2(1,i,j)
           v2(i,j) = ua(i,j)*ec1(2,i,j) + va(i,j)*ec2(2,i,j)
           v3(i,j) = ua(i,j)*ec1(3,i,j) + va(i,j)*ec2(3,i,j)
        enddo
     enddo
     endif
     if ( (je+1)==npy ) then
     do j=npy-1,npy
        do i=is-2,ie+2
           v1(i,j) = ua(i,j)*ec1(1,i,j) + va(i,j)*ec2(1,i,j)
           v2(i,j) = ua(i,j)*ec1(2,i,j) + va(i,j)*ec2(2,i,j)
           v3(i,j) = ua(i,j)*ec1(3,i,j) + va(i,j)*ec2(3,i,j)
        enddo
     enddo
     endif

     im2 = (npx-1)/2;   jm2 = (npy-1)/2
! A -> C (across face averaging taking place here):
! Xdir
     call fill3_4corners(v1, v2, v3, 1)

     if ( is==1 ) then
        i=1
        do j=js-2,je+2
           v1d(1,j) = v1(0,j) + v1(1,j)
           v1d(2,j) = v2(0,j) + v2(1,j)
           v1d(3,j) = v3(0,j) + v3(1,j)
        enddo
        do j=js-1,je+1
        if ( j==0 .or. (j>jm2 .and. j<npy) ) then
           do k=1,3
              w3(k,j) = edge_vect_w(j)*v1d(k,j-1)+(1.-edge_vect_w(j))*v1d(k,j)
!             w3(k,j) = v1d(k,j)
           enddo
        else
           do k=1,3
              w3(k,j) = edge_vect_w(j)*v1d(k,j+1)+(1.-edge_vect_w(j))*v1d(k,j)
!             w3(k,j) = v1d(k,j)
           enddo
        endif
        enddo
! Projection:
        do j=js-1,je+1
           uc(i,j) = 0.5*(w3(1,j)*ew(1,i,j,1)+w3(2,j)*ew(2,i,j,1)+w3(3,j)*ew(3,i,j,1))
           ut(i,j) = uc(i,j)*rsin_u(i,j)
        enddo
     endif

     if ( (ie+1)==npx ) then
        i=npx
        do j=js-2,je+2
           v1d(1,j) = v1(ie,j) + v1(i,j)
           v1d(2,j) = v2(ie,j) + v2(i,j)
           v1d(3,j) = v3(ie,j) + v3(i,j)
        enddo
        do j=js-1,je+1
        if ( j==0 .or. (j>jm2 .and. j<npy) ) then
           do k=1,3
              w3(k,j) = edge_vect_e(j)*v1d(k,j-1)+(1.-edge_vect_e(j))*v1d(k,j)
!             w3(k,j) = v1d(k,j)
           enddo
        else
           do k=1,3
              w3(k,j) = edge_vect_e(j)*v1d(k,j+1)+(1.-edge_vect_e(j))*v1d(k,j)
!             w3(k,j) = v1d(k,j)
           enddo
        endif
        enddo
! Projection:
        do j=js-1,je+1
           uc(i,j) = 0.5*(w3(1,j)*ew(1,i,j,1)+w3(2,j)*ew(2,i,j,1)+w3(3,j)*ew(3,i,j,1))
           ut(i,j) = uc(i,j)*rsin_u(i,j)
        enddo
     endif

! Ydir:
     call fill3_4corners(v1, v2, v3, 2)

     if ( js==1 ) then
        j=1
        do i=is-2,ie+2
           v1d(1,i) = v1(i,0) + v1(i,1)
           v1d(2,i) = v2(i,0) + v2(i,1)
           v1d(3,i) = v3(i,0) + v3(i,1)
        enddo
        do i=is-1,ie+1
        if ( i<=0 .or. (i>im2 .and. i<npx) ) then
           do k=1,3
              w3(k,i) = edge_vect_s(i)*v1d(k,i-1)+ (1.-edge_vect_s(i))*v1d(k,i)
!             w3(k,i) = v1d(k,i)
           enddo
        else
           do k=1,3
              w3(k,i) = edge_vect_s(i)*v1d(k,i+1)+(1.-edge_vect_s(i))*v1d(k,i)
!             w3(k,i) = v1d(k,i)
           enddo
        endif
        enddo
! Projection:
        do i=is-1,ie+1
           vc(i,j) = 0.5*(w3(1,i)*es(1,i,j,2)+w3(2,i)*es(2,i,j,2)+w3(3,i)*es(3,i,j,2))
           vt(i,j) = vc(i,j)*rsin_v(i,j)
        enddo
     endif

     if ( (je+1)==npy ) then
        j=npy
        do i=is-2,ie+2
           v1d(1,i) = v1(i,je) + v1(i,j)
           v1d(2,i) = v2(i,je) + v2(i,j)
           v1d(3,i) = v3(i,je) + v3(i,j)
        enddo
        do i=is-1,ie+1
        if ( i<=0 .or. (i>im2 .and. i<npx) ) then
           do k=1,3
              w3(k,i) = edge_vect_n(i)*v1d(k,i-1)+(1.-edge_vect_n(i))*v1d(k,i)
!             w3(k,i) = v1d(k,i)
           enddo
        else
           do k=1,3
              w3(k,i) = edge_vect_n(i)*v1d(k,i+1)+(1.-edge_vect_n(i))*v1d(k,i)
!             w3(k,i) = v1d(k,i)
           enddo
        endif
        enddo
! Projection:
        do i=is-1,ie+1
           vc(i,j) = 0.5*(w3(1,i)*es(1,i,j,2)+w3(2,i)*es(2,i,j,2)+w3(3,i)*es(3,i,j,2))
           vt(i,j) = vc(i,j)*rsin_v(i,j)
        enddo
     endif

#ifdef VECT_DEV
!  return
!----------
! Final fix
!--------- Use lat-lon component along equator --------------------------
   if ( tile==1 .or. tile==2 ) then
     if ( is==1 ) then
        do j=js,je
           uc(1,j) = a11(0,j)*utmp(0,j)+a12(0,j)*vtmp(0,j) +   &
                     a11(1,j)*utmp(1,j)+a12(1,j)*vtmp(1,j)
           ut(1,j) = uc(1,j)*rsin_u(1,j)
        enddo
     endif
     if ( (ie+1)==npx ) then
        do j=js,je
           uc(npx,j) = a11(ie, j)*utmp(ie, j)+a12(ie, j)*vtmp(ie, j) +   &
                       a11(npx,j)*utmp(npx,j)+a12(npx,j)*vtmp(npx,j)
           ut(npx,j) = uc(npx,j)*rsin_u(npx,j)
        enddo
     endif
   endif

   if ( tile==4 .or. tile==5 ) then
     if ( js==1 ) then
        do i=is,ie
           vc(i,1) = a11(i,0)*utmp(i,0)+a12(i,0)*vtmp(i,0) +   &
                     a11(i,1)*utmp(i,1)+a12(i,1)*vtmp(i,1)
           vt(i,1) = vc(i,1)*rsin_v(i,1)
        enddo
     endif
     if ( (je+1)==npy ) then
        do i=is,ie
           vc(i,npy) = a11(i,je )*utmp(i,je )+a12(i,je )*vtmp(i,je ) +   &
                       a11(i,npy)*utmp(i,npy)+a12(i,npy)*vtmp(i,npy)
           vt(i,npy) = vc(i,npy)*rsin_v(i,npy)
        enddo
     endif
   endif
!  call mpp_update_domains(uc, vc, domain, gridtype=CGRID_NE, complete=.true.)
#endif


 end subroutine d2a2c_vect_v0

 
      
!-------------------------------------------------------------------------------
 subroutine d2a2c_vect_v1( u,v, ua,va, uc,vc, ut,vt )
  real, intent(in) ::  u(isd:ied,jsd:jed+1)
  real, intent(in) ::  v(isd:ied+1,jsd:jed)
  real, intent(out), dimension(isd:ied+1,jsd:jed  ):: uc
  real, intent(out), dimension(isd:ied  ,jsd:jed+1):: vc
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: ua, va, ut, vt
!--------------------------------------------------------------
! Local 
  real, dimension(isd:ied,jsd:jed):: v1, v2, v3
    real wk(is-2:ie+2,js-2:je+2)
    real v1d(3,-3:npx+4), w3(3,-3:npx+4)

    real vw1, vw2, vw3
    real vs1, vs2, vs3
    integer i, j, k, im, jm, im2, jm2
    real utmp, vtmp

     im = npx - 1
     jm = npy - 1

! needs only ut[is-1:ie+2,js-1:je+1], vt[is-1:ie+1,js-1:je+2]

     do j=js-2,je+3
        do i=is-2,ie+2
           vt(i,j) = u(i,j)*dx(i,j)
        enddo
     enddo
     do j=js-2,je+2
        do i=is-2,ie+3
           ut(i,j) = v(i,j)*dy(i,j)
        enddo
     enddo

     do j=js-2,je+2
        do i=is-2,ie+2
! Co-variant to Co-variant "vorticity-conserving" interpolation
           utmp = 0.5*(vt(i,j) + vt(i,j+1)) * rdxa(i,j)
           vtmp = 0.5*(ut(i,j) + ut(i+1,j)) * rdya(i,j)
           ua(i,j) = (utmp-vtmp*cosa_s(i,j)) * rsin2(i,j)
           va(i,j) = (vtmp-utmp*cosa_s(i,j)) * rsin2(i,j)
        enddo
     enddo

! V = ua * e1 + va * e2
     do j=js-2,je+2
        do i=is-2,ie+2
           v1(i,j) = ua(i,j)*ec1(1,i,j) + va(i,j)*ec2(1,i,j)
           v2(i,j) = ua(i,j)*ec1(2,i,j) + va(i,j)*ec2(2,i,j)
           v3(i,j) = ua(i,j)*ec1(3,i,j) + va(i,j)*ec2(3,i,j)
        enddo
     enddo

     im2 = (npx-1)/2;   jm2 = (npy-1)/2

! A -> C (across face averaging taking place here):
! Xdir
     call fill3_4corners(v1, v2, v3, 1)

! A -> C
     do j=js-1,je+1
        do i=is-1,ie+2
           vw1 = v1(i-1,j) + v1(i,j)
           vw2 = v2(i-1,j) + v2(i,j)
           vw3 = v3(i-1,j) + v3(i,j)
           uc(i,j) = 0.5*(vw1*ew(1,i,j,1) + vw2*ew(2,i,j,1) + vw3*ew(3,i,j,1))
              vtmp = 0.5*(vw1*ew(1,i,j,2) + vw2*ew(2,i,j,2) + vw3*ew(3,i,j,2))
           ut(i,j) = (uc(i,j)-vtmp*cosa_u(i,j)) * rsin_u(i,j)
        enddo
     enddo

! Fix the edge:

     if ( is==1 ) then
        i=1
        do j=js-2,je+2
           v1d(1,j) = v1(i-1,j) + v1(i,j)
           v1d(2,j) = v2(i-1,j) + v2(i,j)
           v1d(3,j) = v3(i-1,j) + v3(i,j)
        enddo
        do j=js-1,je+1
        if ( j==0 .or. (j>jm2 .and. j<npy) ) then
           do k=1,3
              w3(k,j) = edge_vect_w(j)*v1d(k,j-1)+(1.-edge_vect_w(j))*v1d(k,j)
           enddo
        else
           do k=1,3
              w3(k,j) = edge_vect_w(j)*v1d(k,j+1)+(1.-edge_vect_w(j))*v1d(k,j)
           enddo
        endif
        enddo
! Projection:
        do j=js-1,je+1
           uc(i,j) = 0.5*(w3(1,j)*ew(1,i,j,1)+w3(2,j)*ew(2,i,j,1)+w3(3,j)*ew(3,i,j,1))
           ut(i,j) = uc(i,j)*rsin_u(i,j)
        enddo
     endif

     if ( (ie+1)==npx ) then
        i=npx
        do j=js-2,je+2
           v1d(1,j) = v1(i-1,j) + v1(i,j)
           v1d(2,j) = v2(i-1,j) + v2(i,j)
           v1d(3,j) = v3(i-1,j) + v3(i,j)
        enddo
        do j=js-1,je+1
        if ( j==0 .or. (j>jm2 .and. j<npy) ) then
           do k=1,3
              w3(k,j) = edge_vect_e(j)*v1d(k,j-1)+(1.-edge_vect_e(j))*v1d(k,j)
           enddo
        else
           do k=1,3
              w3(k,j) = edge_vect_e(j)*v1d(k,j+1)+(1.-edge_vect_e(j))*v1d(k,j)
           enddo
        endif
        enddo
! Projection:
        do j=js-1,je+1
           uc(i,j) = 0.5*(w3(1,j)*ew(1,i,j,1)+w3(2,j)*ew(2,i,j,1)+w3(3,j)*ew(3,i,j,1))
           ut(i,j) = uc(i,j)*rsin_u(i,j)
        enddo
     endif

! Ydir:
     call fill3_4corners(v1, v2, v3, 2)

     do j=js-1,je+2
        do i=is-1,ie+1
           vs1 = v1(i,j-1) + v1(i,j)
           vs2 = v2(i,j-1) + v2(i,j)
           vs3 = v3(i,j-1) + v3(i,j)
           vc(i,j) = 0.5*(vs1*es(1,i,j,2) + vs2*es(2,i,j,2) + vs3*es(3,i,j,2))
              utmp = 0.5*(vs1*es(1,i,j,1) + vs2*es(2,i,j,1) + vs3*es(3,i,j,1))
           vt(i,j) = (vc(i,j)-utmp*cosa_v(i,j)) * rsin_v(i,j)
        enddo
     enddo

! Fix the edge:
     if ( js==1 ) then
        j=1
        do i=is-2,ie+2
           v1d(1,i) = v1(i,j-1) + v1(i,j)
           v1d(2,i) = v2(i,j-1) + v2(i,j)
           v1d(3,i) = v3(i,j-1) + v3(i,j)
        enddo
        do i=is-1,ie+1
        if ( i<=0 .or. (i>im2 .and. i<npx) ) then
           do k=1,3
              w3(k,i) = edge_vect_s(i)*v1d(k,i-1)+ (1.-edge_vect_s(i))*v1d(k,i)
           enddo
        else
           do k=1,3
              w3(k,i) = edge_vect_s(i)*v1d(k,i+1)+(1.-edge_vect_s(i))*v1d(k,i)
           enddo
        endif
        enddo
! Projection:
        do i=is-1,ie+1
           vc(i,j) = 0.5*(w3(1,i)*es(1,i,j,2)+w3(2,i)*es(2,i,j,2)+w3(3,i)*es(3,i,j,2))
           vt(i,j) = vc(i,j)*rsin_v(i,j)
        enddo
     endif

     if ( (je+1)==npy ) then
        j=npy
        do i=is-2,ie+2
           v1d(1,i) = v1(i,j-1) + v1(i,j)
           v1d(2,i) = v2(i,j-1) + v2(i,j)
           v1d(3,i) = v3(i,j-1) + v3(i,j)
        enddo
        do i=is-1,ie+1
        if ( i<=0 .or. (i>im2 .and. i<npx) ) then
           do k=1,3
              w3(k,i) = edge_vect_n(i)*v1d(k,i-1)+(1.-edge_vect_n(i))*v1d(k,i)
           enddo
        else
           do k=1,3
              w3(k,i) = edge_vect_n(i)*v1d(k,i+1)+(1.-edge_vect_n(i))*v1d(k,i)
           enddo
        endif
        enddo
! Projection:
        do i=is-1,ie+1
           vc(i,j) = 0.5*(w3(1,i)*es(1,i,j,2)+w3(2,i)*es(2,i,j,2)+w3(3,i)*es(3,i,j,2))
           vt(i,j) = vc(i,j)*rsin_v(i,j)
        enddo
     endif

 end subroutine d2a2c_vect_v1



 subroutine fill3_4corners(q1, q2, q3, dir)
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
  integer, intent(in):: dir                ! 1: x-dir; 2: y-dir
  real, intent(inout):: q1(isd:ied,jsd:jed)
  real, intent(inout):: q2(isd:ied,jsd:jed)
  real, intent(inout):: q3(isd:ied,jsd:jed)

  select case(dir)
  case(1)
      if ( sw_corner ) then
          q1(-1,0) = q1(0,2); q1(0,0) = q1(0,1); q1(0,-1) = q1(-1,1)
          q2(-1,0) = q2(0,2); q2(0,0) = q2(0,1); q2(0,-1) = q2(-1,1)
          q3(-1,0) = q3(0,2); q3(0,0) = q3(0,1); q3(0,-1) = q3(-1,1)
      endif
      if ( se_corner ) then
          q1(npx+1,0) = q1(npx,2); q1(npx,0) = q1(npx,1); q1(npx,-1) = q1(npx+1,1)
          q2(npx+1,0) = q2(npx,2); q2(npx,0) = q2(npx,1); q2(npx,-1) = q2(npx+1,1)
          q3(npx+1,0) = q3(npx,2); q3(npx,0) = q3(npx,1); q3(npx,-1) = q3(npx+1,1)
      endif
      if ( ne_corner ) then
          q1(npx,npy) = q1(npx,npy-1); q1(npx+1,npy) = q1(npx,npy-2); q1(npx,npy+1) = q1(npx+1,npy-1)
          q2(npx,npy) = q2(npx,npy-1); q2(npx+1,npy) = q2(npx,npy-2); q2(npx,npy+1) = q2(npx+1,npy-1)
          q3(npx,npy) = q3(npx,npy-1); q3(npx+1,npy) = q3(npx,npy-2); q3(npx,npy+1) = q3(npx+1,npy-1)
      endif
      if ( nw_corner ) then
          q1(0,npy) = q1(0,npy-1); q1(-1,npy) = q1(0,npy-2); q1(0,npy+1) = q1(-1,npy-1)
          q2(0,npy) = q2(0,npy-1); q2(-1,npy) = q2(0,npy-2); q2(0,npy+1) = q2(-1,npy-1)
          q3(0,npy) = q3(0,npy-1); q3(-1,npy) = q3(0,npy-2); q3(0,npy+1) = q3(-1,npy-1)
      endif

  case(2)
      if ( sw_corner ) then
          q1(0,0) = q1(1,0); q1(0,-1) = q1(2,0); q1(-1,0) = q1(1,-1)
          q2(0,0) = q2(1,0); q2(0,-1) = q2(2,0); q2(-1,0) = q2(1,-1)
          q3(0,0) = q3(1,0); q3(0,-1) = q3(2,0); q3(-1,0) = q3(1,-1)
      endif
      if ( se_corner ) then
          q1(npx,0) = q1(npx-1,0); q1(npx,-1) = q1(npx-2,0); q1(npx+1,0) = q1(npx-1,-1)
          q2(npx,0) = q2(npx-1,0); q2(npx,-1) = q2(npx-2,0); q2(npx+1,0) = q2(npx-1,-1)
          q3(npx,0) = q3(npx-1,0); q3(npx,-1) = q3(npx-2,0); q3(npx+1,0) = q3(npx-1,-1)
      endif
      if ( ne_corner ) then
          q1(npx,npy) = q1(npx-1,npy); q1(npx,npy+1) = q1(npx-2,npy); q1(npx+1,npy) = q1(npx-1,npy+1)
          q2(npx,npy) = q2(npx-1,npy); q2(npx,npy+1) = q2(npx-2,npy); q2(npx+1,npy) = q2(npx-1,npy+1)
          q3(npx,npy) = q3(npx-1,npy); q3(npx,npy+1) = q3(npx-2,npy); q3(npx+1,npy) = q3(npx-1,npy+1)
      endif
      if ( nw_corner ) then
          q1(0,npy) = q1(1,npy); q1(0,npy+1) = q1(2,npy); q1(-1,npy) = q1(1,npy+1)
          q2(0,npy) = q2(1,npy); q2(0,npy+1) = q2(2,npy); q2(-1,npy) = q2(1,npy+1)
          q3(0,npy) = q3(1,npy); q3(0,npy+1) = q3(2,npy); q3(-1,npy) = q3(1,npy+1)
      endif

  end select
 end subroutine fill3_4corners

 subroutine fill2_4corners(q1, q2, dir)
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
  integer, intent(in):: dir                ! 1: x-dir; 2: y-dir
  real, intent(inout):: q1(isd:ied,jsd:jed)
  real, intent(inout):: q2(isd:ied,jsd:jed)

  select case(dir)
  case(1)
      if ( sw_corner ) then
          q1(-1,0) = q1(0,2);    q1(0,0) = q1(0,1)
          q2(-1,0) = q2(0,2);    q2(0,0) = q2(0,1)
      endif
      if ( se_corner ) then
          q1(npx+1,0) = q1(npx,2); q1(npx,0) = q1(npx,1)
          q2(npx+1,0) = q2(npx,2); q2(npx,0) = q2(npx,1)
      endif
      if ( nw_corner ) then
          q1(0,npy) = q1(0,npy-1); q1(-1,npy) = q1(0,npy-2)
          q2(0,npy) = q2(0,npy-1); q2(-1,npy) = q2(0,npy-2)
      endif
      if ( ne_corner ) then
          q1(npx,npy) = q1(npx,npy-1); q1(npx+1,npy) = q1(npx,npy-2)
          q2(npx,npy) = q2(npx,npy-1); q2(npx+1,npy) = q2(npx,npy-2)
      endif

  case(2)
      if ( sw_corner ) then
          q1(0,0) = q1(1,0); q1(0,-1) = q1(2,0)
          q2(0,0) = q2(1,0); q2(0,-1) = q2(2,0)
      endif
      if ( se_corner ) then
          q1(npx,0) = q1(npx-1,0); q1(npx,-1) = q1(npx-2,0)
          q2(npx,0) = q2(npx-1,0); q2(npx,-1) = q2(npx-2,0)
      endif
      if ( nw_corner ) then
          q1(0,npy) = q1(1,npy); q1(0,npy+1) = q1(2,npy)
          q2(0,npy) = q2(1,npy); q2(0,npy+1) = q2(2,npy)
      endif
      if ( ne_corner ) then
          q1(npx,npy) = q1(npx-1,npy); q1(npx,npy+1) = q1(npx-2,npy)
          q2(npx,npy) = q2(npx-1,npy); q2(npx,npy+1) = q2(npx-2,npy)
      endif

  end select

 end subroutine fill2_4corners

 subroutine fill_4corners(q, dir)
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
  integer, intent(in):: dir                ! 1: x-dir; 2: y-dir
  real, intent(inout):: q(isd:ied,jsd:jed)

  select case(dir)
  case(1)
      if ( sw_corner ) then
          q(-1,0) = q(0,2)
          q( 0,0) = q(0,1)
      endif
      if ( se_corner ) then
          q(npx+1,0) = q(npx,2)
          q(npx,  0) = q(npx,1)
      endif
      if ( nw_corner ) then
          q( 0,npy) = q(0,npy-1)
          q(-1,npy) = q(0,npy-2)
      endif
      if ( ne_corner ) then
          q(npx,  npy) = q(npx,npy-1)
          q(npx+1,npy) = q(npx,npy-2)
      endif

  case(2)
      if ( sw_corner ) then
          q(0, 0) = q(1,0)
          q(0,-1) = q(2,0)
      endif
      if ( se_corner ) then
          q(npx, 0) = q(npx-1,0)
          q(npx,-1) = q(npx-2,0)
      endif
      if ( nw_corner ) then
          q(0,npy  ) = q(1,npy)
          q(0,npy+1) = q(2,npy)
      endif
      if ( ne_corner ) then
          q(npx,npy  ) = q(npx-1,npy)
          q(npx,npy+1) = q(npx-2,npy)
      endif

  end select

 end subroutine fill_4corners

 end module sw_core

