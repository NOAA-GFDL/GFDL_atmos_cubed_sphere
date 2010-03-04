 module fv_grid_utils_mod
 
#include <fms_platform.h>

 use mpp_mod,         only: FATAL, mpp_error
 use mpp_domains_mod, only: mpp_update_domains, DGRID_NE, mpp_global_sum,   &
                            BITWISE_EXACT_SUM
 use mpp_parameter_mod, only: AGRID_PARAM=>AGRID, CGRID_NE_PARAM=>CGRID_NE, & 
                              CORNER, SCALAR_PAIR

 use external_sst_mod, only: i_sst, j_sst, sst_ncep, sst_anom
 use fv_arrays_mod,   only: fv_atmos_type
 use fv_eta_mod,      only: set_eta
 use fv_mp_mod,       only: domain, ng, is,js,ie,je, isd,jsd,ied,jed, gid,  &
                            mp_reduce_sum, mp_reduce_min, mp_reduce_max
 use fv_timing_mod,   only: timing_on, timing_off

 implicit none
 private
#ifdef NO_QUAD_PRECISION
! 64-bit precision (kind=8)
 integer, parameter:: f_p = selected_real_kind(15)
#else
! Higher precision (kind=16) for grid geometrical factors:
 integer, parameter:: f_p = selected_real_kind(20)
#endif
 real, parameter::  big_number=1.E35
 real, parameter:: tiny_number=1.E-35

! For computing mismatch for variable grid zise:
 real, allocatable :: cx1(:,:), cx2(:,:)    !
 real, allocatable :: cy1(:,:), cy2(:,:)    !

! Scalars:
 real, allocatable :: edge_s(:)
 real, allocatable :: edge_n(:)
 real, allocatable :: edge_w(:)
 real, allocatable :: edge_e(:)
! Vector:
 real, allocatable :: edge_vect_s(:)
 real, allocatable :: edge_vect_n(:)
 real, allocatable :: edge_vect_w(:)
 real, allocatable :: edge_vect_e(:)
! scalar:
 real, allocatable :: ex_s(:)
 real, allocatable :: ex_n(:)
 real, allocatable :: ex_w(:)
 real, allocatable :: ex_e(:)
! Vandermonde Matrix:
 real, allocatable :: van2(:,:,:)
! divergence Damping:
 real, allocatable :: divg_u(:,:), divg_v(:,:)    !
! Cubed_2_latlon:
 real, allocatable :: a11(:,:)
 real, allocatable :: a12(:,:)
 real, allocatable :: a21(:,:)
 real, allocatable :: a22(:,:)
! latlon_2_cubed:
 real, allocatable :: z11(:,:)
 real, allocatable :: z12(:,:)
 real, allocatable :: z21(:,:)
 real, allocatable :: z22(:,:)

 real:: global_area, da_min, da_max, da_min_c, da_max_c
 logical:: g_sum_initialized
 logical:: Gnomonic_grid
 logical:: sw_corner, se_corner, ne_corner, nw_corner 
 real, allocatable :: cosa_u(:,:)
 real, allocatable :: cosa_v(:,:)
 real, allocatable :: cosa_s(:,:)
 real, allocatable :: sina_s(:,:)
 real, allocatable :: sina_u(:,:)
 real, allocatable :: sina_v(:,:)
 real, allocatable :: rsin_u(:,:)
 real, allocatable :: rsin_v(:,:)
 real, allocatable ::  rsina(:,:)
 real, allocatable ::  rsin2(:,:)
 real, allocatable :: ee1(:,:,:)
 real, allocatable :: ee2(:,:,:)
 real, allocatable :: ec1(:,:,:)
 real, allocatable :: ec2(:,:,:)
 real, allocatable :: ew(:,:,:,:)
 real, allocatable :: es(:,:,:,:)

! Unit Normal vectors at cell edges:
 real, allocatable :: en1(:,:,:)
 real, allocatable :: en2(:,:,:)

! Extended Cubed cross-edge winds
 real, allocatable :: eww(:,:)
 real, allocatable :: ess(:,:)

! Unit vectors for lat-lon grid
 real, allocatable :: vlon(:,:,:), vlat(:,:,:)
 real, allocatable :: fC(:,:), f0(:,:)
 real :: deglat=15.

 real, parameter:: ptop_min=1.E-8
 real    :: ptop
 integer :: ks
 integer :: g_type, npxx, npyy
 integer :: c2l_ord

 public ptop, ks, ptop_min, fC, f0, deglat, big_number, ew, es, eww, ess, ec1, ec2
 public sina_u, sina_v, cosa_u, cosa_v, cosa_s, sina_s, rsin_u, rsin_v, rsina, rsin2
 public project_sphere_v, latlon2xyz,  gnomonic_grids, global_area,         &
        sw_corner, se_corner, ne_corner, nw_corner, global_mx,              &
        da_min, da_min_c, edge_s, edge_n, edge_w, edge_e,   &
        edge_vect_s,edge_vect_n,edge_vect_w,edge_vect_e, unit_vect_latlon,  &
        cubed_to_latlon, c2l_ord2, g_sum, global_qsum, great_circle_dist,  &
        v_prod, en1, en2, ex_w, ex_e, ex_s, ex_n, vlon, vlat, ee1, ee2, &
        cx1, cx2, cy1, cy2, Gnomonic_grid, van2, divg_u, divg_v
 public mid_pt_sphere,  mid_pt_cart, vect_cross, grid_utils_init, grid_utils_end, &
        spherical_angle, cell_center2, get_area, inner_prod, fill_ghost,    &
        make_eta_level, expand_cell, cart_to_latlon, intp_great_circle, normalize_vect
 public z11, z12, z21, z22

 contains

   subroutine grid_utils_init(Atm, npx, npy, npz, grid, agrid, area, area_c,  &
                              cosa, sina, dx, dy, dxa, dya, dxc, dyc, non_ortho,   &
                              uniform_ppm, grid_type, c2l_order)
! Initialize 2D memory and geometrical factors
      type(fv_atmos_type), intent(inout) :: Atm
      logical, intent(in):: non_ortho
      integer, intent(in):: npx, npy, npz
      integer, intent(in):: grid_type, c2l_order
      real, intent(in)::  grid(isd:ied+1,jsd:jed+1,2)
      real, intent(in):: agrid(isd:ied  ,jsd:jed  ,2)
      real, intent(in):: area(isd:ied,jsd:jed)
      real, intent(in):: area_c(isd:ied+1,jsd:jed+1)
      real, intent(in)::  dx(isd:ied  ,jsd:jed+1)
      real, intent(in)::  dy(isd:ied+1,jsd:jed  )
      real, intent(inout):: dxa(isd:ied  ,jsd:jed  )
      real, intent(inout):: dya(isd:ied  ,jsd:jed  )
      real, intent(inout):: dxc(isd:ied+1,jsd:jed  )
      real, intent(inout):: dyc(isd:ied  ,jsd:jed+1)

      real, intent(inout):: cosa(isd:ied+1,jsd:jed+1)
      real, intent(inout):: sina(isd:ied+1,jsd:jed+1)
      logical, intent(IN) :: uniform_ppm
!
      real grid3(3,isd:ied+1,jsd:jed+1)
      real p1(3), p2(3), p3(3), pp(3)
      real sin2, tmp1, tmp2
      integer i, j, k, n

      npxx = npx;  npyy = npy

      g_sum_initialized = .false.

      allocate ( Atm%ak(npz+1) )
      allocate ( Atm%bk(npz+1) )

      if ( npz == 1 ) then
           Atm%ak(1) = 0.
           Atm%ak(2) = 0.
           Atm%bk(1) = 0.
           Atm%bk(2) = 1.
           ptop      = 0.
           Atm%ks    = 0
      else
! Initialize (ak,bk) for cold start; overwritten with restart file
           call set_eta(npz, ks, ptop, Atm%ak, Atm%bk)
           Atm%ks = ks
#ifdef PRINT_GRID
           if ( gid==0 ) then
                write(*,*) 'Grid_init', ks, ptop
              do k=1,npz
                 write(*,*) k, atm%ak(k), atm%bk(k)
              enddo
           endif
#endif
      endif

! NCEP analysis available from amip-Interp (allocate if needed)
      if (.not. allocated(sst_ncep)) allocate (sst_ncep(i_sst,j_sst))
      if (.not. allocated(sst_anom)) allocate (sst_anom(i_sst,j_sst))

! Coriolis parameters:
      allocate ( f0(isd:ied  ,jsd:jed  ) )
      allocate ( fC(isd:ied+1,jsd:jed+1) )

! Corner unit vectors:
      allocate( ee1(3,isd:ied+1,jsd:jed+1) )
      allocate( ee2(3,isd:ied+1,jsd:jed+1) )

! Center unit vectors:
      allocate( ec1(3,isd:ied,jsd:jed) )
      allocate( ec2(3,isd:ied,jsd:jed) )

! Edge unit vectors:
      allocate( ew(3,isd:ied+1,jsd:jed,  2) )
      allocate( es(3,isd:ied  ,jsd:jed+1,2) )

! Edge unit "Normal" vectors: (for omega computation)
      allocate( en1(3,is:ie,  js:je+1) )   ! E-W edges
      allocate( en2(3,is:ie+1,js:je  ) )   ! N-S egdes
 
      allocate ( cosa_u(isd:ied+1,jsd:jed) )
      allocate ( sina_u(isd:ied+1,jsd:jed) )
      allocate ( rsin_u(isd:ied+1,jsd:jed) )

      allocate ( cosa_v(isd:ied,jsd:jed+1) )
      allocate ( sina_v(isd:ied,jsd:jed+1) )
      allocate ( rsin_v(isd:ied,jsd:jed+1) )

      allocate ( cosa_s(isd:ied,jsd:jed) )    ! cell center
      allocate ( sina_s(isd:ied,jsd:jed) )    ! cell center

      allocate (  rsina(is:ie+1,js:je+1) )    ! cell corners
      allocate (  rsin2(isd:ied,jsd:jed) )    ! cell center

      allocate( eww(3,4) )
      allocate( ess(3,4) )

! For diveregnce damping:
      allocate (  divg_u(isd:ied,  jsd:jed+1) )
      allocate (  divg_v(isd:ied+1,jsd:jed) )

      sw_corner = .false.
      se_corner = .false.
      ne_corner = .false.
      nw_corner = .false.

      if (grid_type < 3) then
         if (       is==1 .and.  js==1 )      sw_corner = .true.
         if ( (ie+1)==npx .and.  js==1 )      se_corner = .true.
         if ( (ie+1)==npx .and. (je+1)==npy ) ne_corner = .true.
         if (       is==1 .and. (je+1)==npy ) nw_corner = .true.
      endif

! For variable grid
   if ( .not. uniform_ppm ) then
      allocate ( cx1(isd:ied,jsd:jed) )
      allocate ( cx2(isd:ied,jsd:jed) )
      allocate ( cy1(isd:ied,jsd:jed) )
      allocate ( cy2(isd:ied,jsd:jed) )

     do j=jsd,jed
        do i=is-2,ie+2
               tmp1 = dxa(i,j)/(dxa(i-1,j) + dxa(i,j) + dxa(i+1,j))
           cx1(i,j) = tmp1*(dxa(i+1,j)+0.5*dxa(i,j))/(dxa(i-1,j)+dxa(i,j)) 
           cx2(i,j) = tmp1*(dxa(i-1,j)+0.5*dxa(i,j))/(dxa(i,j)+dxa(i+1,j)) 
        enddo
     enddo

     do j=js-2,je+2
        do i=isd,ied
               tmp2 = dya(i,j)/(dya(i,j-1) + dya(i,j) + dya(i,j+1)) 
           cy1(i,j) = tmp2*(dya(i,j+1)+0.5*dya(i,j))/(dya(i,j-1)+dya(i,j))
           cy2(i,j) = tmp2*(dya(i,j-1)+0.5*dya(i,j))/(dya(i,j)+dya(i,j+1)) 
        enddo
     enddo
  endif


  if (grid_type < 3) then
     do j=jsd,jed+1
        do i=isd,ied+1
           call latlon2xyz(grid(i,j,1:2), grid3(1,i,j))
        enddo
     enddo

     call get_center_vect( npx, npy, grid3, ec1, ec2 )

! Fill arbitrary values in the non-existing corner regions:
     do k=1,3
        call fill_ghost(ec1(k,:,:), npx, npy, big_number)
        call fill_ghost(ec2(k,:,:), npx, npy, big_number)
     enddo

     do j=jsd,jed
        do i=isd+1,ied
        if ( (i<1   .and. j<1  ) .or. (i>npx .and. j<1  ) .or.  &
             (i>npx .and. j>(npy-1)) .or. (i<1   .and. j>(npy-1)) ) then
!            (i>npx .and. j>npy) .or. (i<1   .and. j>npy) ) then
             ew(1:3,i,j,1:2) = 0.
        else
           call mid_pt_cart( grid(i,j,1:2), grid(i,j+1,1:2), pp)
           if (i==1) then
              call latlon2xyz( agrid(i,j,1:2), p1)
              call vect_cross(p2, pp, p1)
           elseif(i==npx) then
              call latlon2xyz( agrid(i-1,j,1:2), p1)
              call vect_cross(p2, p1, pp)
           else
              call latlon2xyz( agrid(i-1,j,1:2), p3)
              call latlon2xyz( agrid(i,  j,1:2), p1)
              call vect_cross(p2, p3, p1)
           endif
           call vect_cross(ew(1,i,j,1), p2, pp)
           call normalize_vect(ew(1,i,j,1))
!---
           call vect_cross(p1, grid3(1,i,j), grid3(1,i,j+1))
           call vect_cross(ew(1,i,j,2), p1, pp)
           call normalize_vect(ew(1,i,j,2))
        endif
        enddo
     enddo

     do j=jsd+1,jed
        do i=isd,ied
        if ( (i<1   .and. j<1  ) .or. (i>(npx-1) .and. j<1  ) .or.  &
             (i>(npx-1) .and. j>npy) .or. (i<1   .and. j>npy) ) then
             es(1:3,i,j,1:2) = 0.
        else
           call mid_pt_cart(grid(i,j,1:2), grid(i+1,j,1:2), pp)
           if (j==1) then
              call latlon2xyz( agrid(i,j,1:2), p1)
              call vect_cross(p2, pp, p1)
           elseif (j==npy) then
              call latlon2xyz( agrid(i,j-1,1:2), p1)
              call vect_cross(p2, p1, pp)
           else 
              call latlon2xyz( agrid(i,j  ,1:2), p1)
              call latlon2xyz( agrid(i,j-1,1:2), p3)
              call vect_cross(p2, p3, p1)
           endif
           call vect_cross(es(1,i,j,2), p2, pp)
           call normalize_vect(es(1,i,j,2))
!---
           call vect_cross(p3, grid3(1,i,j), grid3(1,i+1,j))
           call vect_cross(es(1,i,j,1), p3, pp)
           call normalize_vect(es(1,i,j,1))
        endif
        enddo
     enddo
  else
     ec1(1,:,:)=1.
     ec1(2,:,:)=0.
     ec1(3,:,:)=0.

     ec2(1,:,:)=0.
     ec2(2,:,:)=1.
     ec2(3,:,:)=0.

     ew(1,:,:,1)=1.
     ew(2,:,:,1)=0.
     ew(3,:,:,1)=0.
                                   
     ew(1,:,:,2)=0.
     ew(2,:,:,2)=1.
     ew(3,:,:,2)=0.

     es(1,:,:,1)=1.
     es(2,:,:,1)=0.
     es(3,:,:,1)=0.
                                   
     es(1,:,:,2)=0.
     es(2,:,:,2)=1.
     es(3,:,:,2)=0.
  endif

     if ( non_ortho ) then
           cosa_u = big_number
           cosa_v = big_number
           cosa_s = big_number
           sina_s = big_number
           sina_u = big_number
           sina_v = big_number
           rsin_u = big_number
           rsin_v = big_number
           rsina  = big_number
           rsin2  = big_number
#ifndef GLOBAL_TRIG
          cosa = big_number
          sina = big_number
!       if ( Gnomonic_grid ) then
! The following works well ONLY with the Gnomonic grid (type 0)
!           do j=js,je+1
!              do i=is,ie+1
!                 tmp1 = cos_angle(grid3(1,i,j), grid3(1,i+1,j), grid3(1,i,j+1))
!                 tmp2 = cos_angle(grid3(1,i,j), grid3(1,i-1,j), grid3(1,i,j-1))
!                 cosa(i,j) = 0.5*(tmp1+tmp2)
!                 sina(i,j) = sqrt( max(0.,1. -cosa(i,j)**2) )
!              enddo
!          enddo
!       endif

        do j=js,je+1
           do i=is,ie+1
! unit vect in X-dir: ee1
              if (i==1) then
                  call vect_cross(pp, grid3(1,i,  j), grid3(1,i+1,j))
              elseif(i==npx) then
                  call vect_cross(pp, grid3(1,i-1,j), grid3(1,i,  j))
              else
                  call vect_cross(pp, grid3(1,i-1,j), grid3(1,i+1,j))
              endif
              call vect_cross(ee1(1,i,j), pp, grid3(1,i,j))
              call normalize_vect( ee1(1,i,j) )

! unit vect in Y-dir: ee2
              if (j==1) then
                  call vect_cross(pp, grid3(1,i,j  ), grid3(1,i,j+1))
              elseif(j==npy) then
                  call vect_cross(pp, grid3(1,i,j-1), grid3(1,i,j  ))
              else
                  call vect_cross(pp, grid3(1,i,j-1), grid3(1,i,j+1))
              endif
              call vect_cross(ee2(1,i,j), pp, grid3(1,i,j))
              call normalize_vect( ee2(1,i,j) )

              tmp1 = inner_prod(ee1(1,i,j), ee2(1,i,j))
              cosa(i,j) = sign(min(1., abs(tmp1)), tmp1)
              sina(i,j) = sqrt(max(0.,1. -cosa(i,j)**2))
           enddo
        enddo
#endif

! call mpp_update_domains(cosa, domain, position=CORNER)
! The above does not work because cosa at edges should have two values (left and right)

      do j=jsd,jed
         do i=isd+1,ied
                   tmp1 = inner_prod(ew(1,i,j,1), ew(1,i,j,2))
            cosa_u(i,j) = sign( min(1., abs(tmp1)), tmp1 )
            sin2 = 1. - cosa_u(i,j)**2
            sin2 = min(1., sin2)
            sin2 = max(tiny_number, sin2)  ! sin(alpha)**2 >= 0.75
 
            sina_u(i,j) = sqrt( sin2 )
            rsin_u(i,j) =  1. / sin2
         enddo
      enddo

      do j=jsd+1,jed
         do i=isd,ied
                   tmp1 = inner_prod(es(1,i,j,1), es(1,i,j,2))
            cosa_v(i,j) = sign( min(1., abs(tmp1)), tmp1 )
            sin2 = 1. - cosa_v(i,j)**2
            sin2 = min(1., sin2)
            sin2 = max(tiny_number, sin2)
            sina_v(i,j) = sqrt( sin2 )
            rsin_v(i,j) =  1. / sin2
         enddo
      enddo
     
      do j=jsd,jed
         do i=isd,ied
                  tmp1  = inner_prod(ec1(1,i,j), ec2(1,i,j))
            cosa_s(i,j) = sign(min(1., abs(tmp1)), tmp1 )
            sin2 = 1. - cosa_s(i,j)**2
            sin2 = min(1., sin2)
            sin2 = max(tiny_number, sin2)
            sina_s(i,j) = min(1., sqrt(sin2))
            rsin2(i,j) = 1. / sin2
         enddo
      enddo
! Force the model to fail if incorrect corner values are to be used:
     call fill_ghost(cosa_s, npx, npy,  big_number)
     call fill_ghost(sina_s, npx, npy, tiny_number)

!------------------------------------
! Set special sin values at edges:
!------------------------------------
      do j=js,je+1
         do i=is,ie+1
            if ( i==1 .or. i==npx .or. j==1 .or. j==npy ) then
                 rsina(i,j) = 1. / sina(i,j)
            else
                 rsina(i,j) = 1. / sina(i,j)**2
            endif
         enddo
      enddo

      do j=jsd,jed
         do i=isd+1,ied
            if ( i==1 .or. i==npx ) then
                 rsin_u(i,j) = 1. / sina_u(i,j)
            endif
         enddo
      enddo

      do j=jsd+1,jed
         do i=isd,ied
            if ( j==1 .or. j==npy ) then
                 rsin_v(i,j) = 1. / sina_v(i,j)
            endif
         enddo
      enddo

   else
           sina = 1.
           cosa = 0.
           rsina  = 1.
           rsin2  = 1.
           sina_u = 1.
           sina_v = 1.
           cosa_u = 0.        
           cosa_v = 0.        
           cosa_s = 0.        
           sina_s = 1.        
           rsin_u = 1.
           rsin_v = 1.
   endif

   if ( grid_type < 3 ) then

#ifdef USE_NORM_VECT
!-------------------------------------------------------------
! Make normal vect at face edges after consines are computed:
!-------------------------------------------------------------
! for old d2a2c_vect routines
     do j=js-1,je+1
        if ( is==1 ) then
             i=1
             call vect_cross(ew(1,i,j,1), grid3(1,i,j+1), grid3(1,i,j)) 
             call normalize_vect( ew(1,i,j,1) )
        endif
        if ( (ie+1)==npx ) then
             i=npx
             call vect_cross(ew(1,i,j,1), grid3(1,i,j+1), grid3(1,i,j)) 
             call normalize_vect( ew(1,i,j,1) )
        endif
     enddo

     if ( js==1 ) then
        j=1
        do i=is-1,ie+1
             call vect_cross(es(1,i,j,2), grid3(1,i,j),grid3(1,i+1,j)) 
             call normalize_vect( es(1,i,j,2) )
        enddo
     endif
     if ( (je+1)==npy ) then
        j=npy
        do i=is-1,ie+1
             call vect_cross(es(1,i,j,2), grid3(1,i,j),grid3(1,i+1,j)) 
             call normalize_vect( es(1,i,j,2) )
        enddo
     endif
#endif

! For omega computation:
! Unit vectors:
     do j=js,je+1
        do i=is,ie
           call vect_cross(en1(1,i,j), grid3(1,i,j), grid3(1,i+1,j))
           call normalize_vect( en1(1,i,j) )
        enddo
     enddo
     do j=js,je
        do i=is,ie+1
           call vect_cross(en2(1,i,j), grid3(1,i,j+1), grid3(1,i,j)) 
           call normalize_vect( en2(1,i,j) )
        enddo
     enddo
!-------------------------------------------------------------
! Make unit vectors for the coordinate extension:
!-------------------------------------------------------------
#ifdef EXTEND_VG
     if ( sw_corner ) then
        do k=1,3
           ess(k,1) = grid3(k,1,1) - grid3(k,0,2)
        enddo
        call normalize_vect( ess(1,1) )
        do k=1,3
           eww(k,1) = grid3(k,1,1) - grid3(k,2,0)
        enddo
        call normalize_vect( eww(1,1) )
     endif
     if ( se_corner ) then
        do k=1,3
           ess(k,2) = grid3(k,npx+1,2) - grid3(k,npx,1)
        enddo
        call normalize_vect( ess(1,2) )
        do k=1,3
           eww(k,2) = grid3(k,npx,1) - grid3(k,npx-1,0)
        enddo
        call normalize_vect( eww(1,2) )
     endif
     if ( ne_corner ) then
        do k=1,3
           ess(k,3) = grid3(k,npx+1,npy-1) - grid3(k,npx,npy)
        enddo
        call normalize_vect( ess(1,3) )
        do k=1,3
           eww(k,3) = grid3(k,npx-1,npy+1) - grid3(k,npx,npy)
        enddo
        call normalize_vect( eww(1,3) )
     endif
     if ( nw_corner ) then
        do k=1,3
           ess(k,4) = grid3(k,1,npy) - grid3(k,0,npy-1)
        enddo
        call normalize_vect( ess(1,4) )
        do k=1,3
           eww(k,4) = grid3(k,2,npy+1) - grid3(k,1,npy)
        enddo
        call normalize_vect( eww(1,4) )
     endif
#endif
  endif
 
  do j=jsd,jed+1
     do i=isd,ied
        divg_u(i,j) = sina_v(i,j)*dyc(i,j)/dx(i,j)
     enddo
  enddo
  do j=jsd,jed
     do i=isd,ied+1
        divg_v(i,j) = sina_u(i,j)*dxc(i,j)/dy(i,j)
     enddo
  enddo

! Initialize cubed_sphere to lat-lon transformation:
     call init_cubed_to_latlon( agrid, grid_type, c2l_order )

     call global_mx(area, ng, da_min, da_max)
     if( gid==0 ) write(6,*) 'da_max/da_min=', da_max/da_min

     call global_mx_c(area_c(is:ie,js:je), is, ie, js, je, da_min_c, da_max_c)

     if( gid==0 ) write(6,*) 'da_max_c/da_min_c=', da_max_c/da_min_c

!------------------------------------------------
! Initialization for interpolation at face edges
!------------------------------------------------
! A->B scalar:
     if (grid_type < 3) then
        call mpp_update_domains(divg_v, divg_u, domain, flags=SCALAR_PAIR,      &
                                gridtype=CGRID_NE_PARAM, complete=.true.)
        call edge_factors (non_ortho, grid, agrid, npx, npy)
        call efactor_a2c_v(non_ortho, grid, agrid, npx, npy)
!       call extend_cube_s(non_ortho, grid, agrid, npx, npy, .false.)
!       call van2d_init(grid, agrid, npx, npy)
     else
        allocate ( edge_s(npx) )
        allocate ( edge_n(npx) )
        allocate ( edge_w(npy) )
        allocate ( edge_e(npy) )

        allocate ( edge_vect_s(isd:ied) )
        allocate ( edge_vect_n(isd:ied) )
        allocate ( edge_vect_w(jsd:jed) )
        allocate ( edge_vect_e(jsd:jed) )

        allocate ( ex_s(npx) )
        allocate ( ex_n(npx) )
        allocate ( ex_w(npy) )
        allocate ( ex_e(npy) )

        edge_s = big_number
        edge_n = big_number
        edge_w = big_number
        edge_e = big_number

        edge_vect_s = big_number
        edge_vect_n = big_number
        edge_vect_w = big_number
        edge_vect_e = big_number

        ex_s(npx) = big_number
        ex_n(npx) = big_number
        ex_w(npy) = big_number
        ex_e(npy) = big_number

     endif

  end subroutine grid_utils_init

 
  subroutine grid_utils_end(uniform_ppm)
  logical, intent(IN) :: uniform_ppm
 
! deallocate sst_ncep (if allocated)
      if (allocated(sst_ncep)) deallocate( sst_ncep )
      if (allocated(sst_anom)) deallocate( sst_anom )

  if ( .not. uniform_ppm ) then
      deallocate( cx1 )
      deallocate( cx2 )
      deallocate( cy1 )
      deallocate( cy2 )
  endif

      deallocate( cosa_u )
      deallocate( cosa_v )
      deallocate( cosa_s )
      deallocate( sina_s )
      deallocate( sina_u )
      deallocate( sina_v )

      deallocate( rsin_u )
      deallocate( rsin_v )
      deallocate( rsina  )
      deallocate( rsin2  )

!     if ( .not. Gnomonic_grid ) then
           deallocate( ee1 )
           deallocate( ee2 )
!     endif

      deallocate( ec1 )
      deallocate( ec2 )
      deallocate( ew )
      deallocate( es )

      deallocate( en1 )
      deallocate( en2 )

      deallocate( eww )
      deallocate( ess )

      deallocate( edge_s )
      deallocate( edge_n )
      deallocate( edge_w )
      deallocate( edge_e )

      deallocate( edge_vect_s )
      deallocate( edge_vect_n )
      deallocate( edge_vect_w )
      deallocate( edge_vect_e )

      if ( allocated(ex_s) ) then
           deallocate( ex_s )
           deallocate( ex_n )
           deallocate( ex_w )
           deallocate( ex_e )
      endif
      if ( allocated(van2) ) deallocate( van2 )

    if ( g_type<4 ) then
      deallocate( a11 )
      deallocate( a12 )
      deallocate( a21 )
      deallocate( a22 )
      deallocate( vlon )
      deallocate( vlat )
      deallocate( z11 )
      deallocate( z12 )
      deallocate( z21 )
      deallocate( z22 )
    endif

    deallocate( divg_u )
    deallocate( divg_v )

  end subroutine grid_utils_end



  real function inner_prod(v1, v2)
       real ,intent(in):: v1(3), v2(3)
       real (f_p) :: vp1(3), vp2(3), prod16
       integer k
      
         do k=1,3
            vp1(k) = v1(k)
           vp2(k) = v2(k)
         enddo
         prod16 = vp1(1)*vp2(1) + vp1(2)*vp2(2) + vp1(3)*vp2(3)
         inner_prod = prod16

  end function inner_prod


 subroutine van2d_init(grid, agrid0, npx, npy)
  integer, intent(in):: npx, npy
  real,    intent(in)::  grid(isd:ied+1,jsd:jed+1,2)
  real,    intent(in):: agrid0(isd:ied ,jsd:jed  ,2)
!
  integer, parameter:: n16 = 16
  real:: agrid(is-2:ie+2 ,js-2:je+2,2)
  real:: a(n16,n16), b(n16,n16), x(n16), y(n16)
  real:: x3, x2, x1, y3, y2, y1, lat, lon, lat0, lon0, sum0
  real:: cos_lat, sin_lat, cos_lat0, sin_lat0, cosc, mfactor
  integer i, j, k, ip, jp

  do j=js-2, je+2
     do i=is-2, ie+2
        agrid(i,j,1) = agrid0(i,j,1)
        agrid(i,j,2) = agrid0(i,j,2)
     enddo
  enddo

  allocate ( van2(n16,is:ie+1,js:je+1) )

  van2 = 0.
  do 2500 j=js, je+1
     do 2500 i=is, ie+1
            lon0 = grid(i,j,1)
            lat0 = grid(i,j,2)
        cos_lat0 = cos(lat0)
        sin_lat0 = sin(lat0)

!--------------
! fill corners:
!--------------
! SW:
        if ( i==1 .and. j==1 ) go to 2000       ! 12-pt matrix
        if ( i==2 .and. j==1 ) then 
!rab             go to 2000
! shift the commom point
             agrid(-1,-1,1:2) = agrid(2,-1,1:2)    ! k=1
             agrid( 0,-1,1:2) = agrid(2, 0,1:2)    ! k=2
             agrid(-1, 0,1:2) = agrid(1,-1,1:2)    ! k=3
             agrid( 0, 0,1:2) = agrid(1, 0,1:2)    ! k=4
        endif
! shift the commom point
        if ( i==1 .and. j==2 ) goto 2000
        if ( i==2 .and. j==2 ) agrid(0,0,1:2) = agrid(4,4,1:2) ! k=1 

! SE:
        if ( i==npx-1 .and. j==1 ) goto 2000
        if ( i==npx   .and. j==1 ) goto 2000    ! 12-pt matrix
        if ( i==npx-1 .and. j==2 ) agrid(npx,0,1:2) = agrid(npx-4,4,1:2) ! k=4
        if ( i==npx   .and. j==2 ) goto 2000

! NE:
        if ( i==npx-1 .and. j==npy-1) agrid(npx,npy,1:2) = agrid(npx-4,npy-4,1:2) ! k=16
        if ( i==npx   .and. j==npy-1) goto 2000
        if ( i==npx-1 .and. j==npy)   goto 2000
        if ( i==npx   .and. j==npy )  goto 2000 ! 12-pt matrix

! NW:
        if ( i==1 .and. j==npy-1 ) goto 2000
        if ( i==2 .and. j==npy-1 ) agrid(0,npy,1:2) = agrid(4,npy-4,1:2) ! k=13
        if ( i==1 .and. j==npy )   goto 2000     ! 12-pt matrix
        if ( i==2 .and. j==npy )   goto 2000

        do k=1,n16
           if    ( k==1 ) then
                               ip = i-2; jp = j-2
           elseif( k==2 ) then
                               ip = i-1; jp = j-2
           elseif( k==3 ) then
                               ip = i;   jp = j-2
           elseif( k==4 ) then
                               ip = i+1; jp = j-2
           elseif( k==5 ) then
                               ip = i-2; jp = j-1
           elseif( k==6 ) then
                               ip = i-1; jp = j-1
           elseif( k==7 ) then
                               ip = i  ; jp = j-1
           elseif( k==8 ) then
                               ip = i+1; jp = j-1
           elseif( k==9 ) then
                               ip = i-2; jp = j
           elseif( k==10 ) then
                               ip = i-1; jp = j
           elseif( k==11 ) then
                               ip = i;   jp = j
           elseif( k==12 ) then
                               ip = i+1; jp = j
           elseif( k==13 ) then
                               ip = i-2; jp = j+1
           elseif( k==14 ) then
                               ip = i-1; jp = j+1
           elseif( k==15 ) then
                               ip = i;   jp = j+1
           elseif( k==16 ) then
                               ip = i+1; jp = j+1
           endif

           lon = agrid(ip,jp,1) 
           lat = agrid(ip,jp,2) 
  
           cos_lat = cos(lat)
           sin_lat = sin(lat)
! Gnomonic projection:
           mfactor = 1. / (sin_lat*sin_lat0 + cos_lat*cos_lat0*cos(lon-lon0))
           x(k) =  cos_lat *sin(lon-lon0)*mfactor
           y(k) = (cos_lat0*sin_lat-sin_lat0*cos_lat*cos(lon-lon0))*mfactor
        enddo

        do k=1,n16
!-------------------------------------
! Full 16x16 "Vandermonde" Matrix
!-------------------------------------
           x1 = x(k)
           x2 = x1*x1
           x3 = x1*x2
           y1 = y(k)
           y2 = y1*y1
           y3 = y1*y2
           a( 1,k) = x3 * y3
           a( 2,k) = x3 * y2
           a( 3,k) = x2 * y3
           a( 4,k) = x2 * y2
           a( 5,k) = x3 * y1
           a( 6,k) = x2 * y1
           a( 7,k) = x1 * y3
           a( 8,k) = x1 * y2
           a( 9,k) = x1 * y1
           a(10,k) = x3
           a(11,k) = x2
           a(12,k) = x1
           a(13,k) = y3
           a(14,k) = y2
           a(15,k) = y1
           a(16,k) = 1.
        enddo

        call invert_matrix(n16, a, b)

        do k=1,n16
           van2(k,i,j) = b(k,n16)
        enddo

        sum0 = 0.
        do k=1,n16
           sum0 = sum0 + b(k,n16)
#ifdef CHECK_VAN2
           if ( k==1 .and. i==3 .and. j==3 ) then
                write(*,*) k,'Van2(3,3):', van2(k,i,j)
!               write(*,*) '          ', lon0, lat0
           endif
#endif
        enddo
        if (abs(sum0-1.)>1.e-10) call mpp_error(FATAL, 'van2_init')
2000 continue
2500 continue


 end subroutine van2d_init

  subroutine van2_init(xs, ys, npx, npy)
  integer, intent(in):: npx, npy
  real,    intent(in), dimension(npx,npy):: xs, ys   ! coner positions
! Local:
  real, dimension(npx,npy):: lon2, lat2
  real::  grid(isd:ied+1,jsd:jed+1,2)
  real:: agrid(is-2:ie+2,js-2:je+2,2)
  integer, parameter:: n16 = 16
  real:: a(n16,n16), b(n16,n16), x(n16), y(n16)
  real:: x3, x2, x1, y3, y2, y1, lat, lon, lat0, lon0, sum0, xk
  real:: cos_lat, sin_lat, cos_lat0, sin_lat0, cosc, mfactor
  real pi
  integer i, j, k, ip, jp

  pi = 4.*atan(1.)

  do j=1,npy
     do i=1,npx
        lat2(i,j) = ys(i,j)
        lon2(i,j) = xs(i,j)
        if ( lon2(i,j) < 0. ) lon2(i,j) = lon2(i,j) + 2.*pi
     enddo
  enddo

  do j=max(1,jsd), min(npy,jed+1)
     do i=max(1,isd), min(npx,ied+1)
        grid(i,j,1) = lon2(i,j)
        grid(i,j,2) = lat2(i,j)
     enddo
  enddo

! agrid = 0.
  do j=max(1,js-2), min(npy-1,je+2)
     do i=max(1,is-2), min(npx-1,ie+2)
        call cell_center2( grid(i,j,  1:2), grid(i+1,j,  1:2),                &
                           grid(i,j+1,1:2), grid(i+1,j+1,1:2), agrid(i,j,1:2) )
     enddo
  enddo

! Fill outer edges
  if ( is==1 ) then 
       do j=max(1,js-2), min(npy-1,je+2)
          call mirror_latlon(lon2(1,1), lat2(1,1), lon2(1,npy), lat2(1,npy),   &
                             agrid(1,j,1), agrid(1,j,2), agrid(0,j,1), agrid(0,j,2))
          call mirror_latlon(lon2(1,1), lat2(1,1), lon2(1,npy), lat2(1,npy),   &
                             agrid(2,j,1), agrid(2,j,2), agrid(-1,j,1), agrid(-1,j,2))
       enddo
  endif
  if ( (ie+1)==npx ) then 
       do j=max(1,js-2), min(npy-1,je+2)
          call mirror_latlon(lon2(npx,1), lat2(npx,1), lon2(npx,npy), lat2(npx,npy),   &
                             agrid(npx-1,j,1), agrid(npx-1,j,2), agrid(npx,j,1), agrid(npx,j,2))
          call mirror_latlon(lon2(npx,1), lat2(npx,1), lon2(npx,npy), lat2(npx,npy),   &
                             agrid(npx-2,j,1), agrid(npx-2,j,2), agrid(npx+1,j,1), agrid(npx+1,j,2))
       enddo
  endif
  if ( js==1 ) then 
       do i=max(1,is-2), min(npx-1,ie+2)
          call mirror_latlon(lon2(1,1), lat2(1,1), lon2(npx,1), lat2(npx,1),   &
                             agrid(i,1,1), agrid(i,1,2), agrid(i,0,1), agrid(i,0,2))
          call mirror_latlon(lon2(1,1), lat2(1,1), lon2(npx,1), lat2(npx,1),   &
                             agrid(i,2,1), agrid(i,2,2), agrid(i,-1,1), agrid(i,-1,2))
       enddo
  endif
  if ( (je+1)==npy ) then 
       do i=max(1,is-2), min(npx-1,ie+2)
          call mirror_latlon(lon2(1,npy), lat2(1,npy), lon2(npx,npy), lat2(npx,npy),   &
                             agrid(i,npy-1,1), agrid(i,npy-1,2), agrid(i,npy,1), agrid(i,npy,2))
          call mirror_latlon(lon2(1,npy), lat2(1,npy), lon2(npx,npy), lat2(npx,npy),   &
                             agrid(i,npy-2,1), agrid(i,npy-2,2), agrid(i,npy+1,1), agrid(i,npy+1,2))
       enddo
  endif

  allocate ( van2(n16,is:ie+1,js:je+1) )

  van2 = 0.
  do 2500 j=js, je+1
     do 2500 i=is, ie+1
            lon0 = grid(i,j,1)
            lat0 = grid(i,j,2)
        cos_lat0 = cos(lat0)
        sin_lat0 = sin(lat0)
!----
! SW:
!----
        if ( i==1 .and. j==1 ) then 
             go to 2000
        endif
        if ( i==2 .and. j==1 ) then 
             agrid(0,-1,1:2) = agrid(-1,2,1:2)
             agrid(0, 0,1:2) = agrid(-1,1,1:2)
        endif
        if ( i==1 .and. j==2 ) then 
             agrid(-1,0,1:2) = agrid(2,-1,1:2)
             agrid( 0,0,1:2) = agrid(1,-1,1:2)
        endif
        if ( i==2 .and. j==2 ) then 
             agrid(0,0,1:2) = agrid(4,4,1:2)   ! add extra point to make it 16
        endif
!----
! SE:
!----
        if ( i==npx   .and. j==1 ) then 
             go to 2000
        endif
        if ( i==npx-1 .and. j==1 ) then 
             agrid(npx,-1,1:2) = agrid(npx+1,2,1:2)
             agrid(npx, 0,1:2) = agrid(npx+1,1,1:2)
        endif
        if ( i==npx-1 .and. j==2 ) then 
             agrid(npx,0,1:2) = agrid(npx-4,4,1:2)
        endif
        if ( i==npx   .and. j==2 ) then 
             agrid(npx+1,0,1:2) = agrid(npx-2,-1,1:2)
             agrid(npx,  0,1:2) = agrid(npx-1,-1,1:2)
        endif
!----
! NE:
!----
        if ( i==npx   .and. j==npy ) then 
             go to 2000
        endif
        if ( i==npx-1 .and. j==npy-1) then 
             agrid(npx,npy,1:2) = agrid(npx-4,npy-4,1:2)
        endif
        if ( i==npx   .and. j==npy-1) then 
             agrid(npx+1,npy,1:2) = agrid(npx-2,npy+1,1:2)
             agrid(npx,  npy,1:2) = agrid(npx-1,npy+1,1:2)
        endif
        if ( i==npx-1 .and. j==npy) then 
             agrid(npx,npy+1,1:2) = agrid(npx+1,npy-2,1:2)
             agrid(npx,npy,  1:2) = agrid(npx+1,npy-1,1:2)
        endif
!----
! NW:
!----
        if ( i==1 .and. j==npy ) then 
             go to 2000
        endif
        if ( i==1 .and. j==npy-1 ) then 
             agrid(-1,npy,1:2) = agrid(2,npy+1,1:2)
             agrid( 0,npy,1:2) = agrid(1,npy+1,1:2)
        endif
        if ( i==2 .and. j==npy-1 ) then 
             agrid(0,npy,1:2) = agrid(4,npy-4,1:2)
        endif
        if ( i==2 .and. j==npy ) then 
             agrid(0,npy+1,1:2) = agrid(-1,npy-2,1:2)
             agrid(0,npy,  1:2) = agrid(-1,npy-1,1:2)
        endif

        do k=1,n16
           if    ( k==1 ) then
                               ip = i-2; jp = j-2
           elseif( k==2 ) then
                               ip = i-1; jp = j-2
           elseif( k==3 ) then
                               ip = i;   jp = j-2
           elseif( k==4 ) then
                               ip = i+1; jp = j-2
           elseif( k==5 ) then
                               ip = i-2; jp = j-1
           elseif( k==6 ) then
                               ip = i-1; jp = j-1
           elseif( k==7 ) then
                               ip = i  ; jp = j-1
           elseif( k==8 ) then
                               ip = i+1; jp = j-1
           elseif( k==9 ) then
                               ip = i-2; jp = j
           elseif( k==10 ) then
                               ip = i-1; jp = j
           elseif( k==11 ) then
                               ip = i;   jp = j
           elseif( k==12 ) then
                               ip = i+1; jp = j
           elseif( k==13 ) then
                               ip = i-2; jp = j+1
           elseif( k==14 ) then
                               ip = i-1; jp = j+1
           elseif( k==15 ) then
                               ip = i;   jp = j+1
           elseif( k==16 ) then
                               ip = i+1; jp = j+1
           endif

           lon = agrid(ip,jp,1) 
           lat = agrid(ip,jp,2) 
  
           cos_lat = cos(lat)
           sin_lat = sin(lat)
! Gnomonic projection:
           mfactor = 1. / (sin_lat*sin_lat0 + cos_lat*cos_lat0*cos(lon-lon0))
           x(k) =  cos_lat *sin(lon-lon0)*mfactor
           y(k) = (cos_lat0*sin_lat - sin_lat0*cos_lat*cos(lon-lon0))*mfactor
#ifdef MIRROR_V
           if ( j==1 .or. j==npy ) then
                  xk = x(k)
                x(k) = y(k)
                y(k) = xk
           endif
#endif

        enddo

        do k=1,n16
!-------------------------------------
! Full 16x16 "Vandermonde" Matrix
!-------------------------------------
           x1 = x(k)
           x2 = x1*x1
           x3 = x1*x2
           y1 = y(k)
           y2 = y1*y1
           y3 = y1*y2
!---------------------
           a( 1,k) = x3 * y3
           a( 2,k) = x3 * y2
           a( 3,k) = x2 * y3
           a( 4,k) = x2 * y2
           a( 5,k) = x3 * y1
           a( 6,k) = x2 * y1
           a( 7,k) = x1 * y3
           a( 8,k) = x1 * y2
           a( 9,k) = x1 * y1
           a(10,k) = x3
           a(11,k) = x2
           a(12,k) = x1
           a(13,k) = y3
           a(14,k) = y2
           a(15,k) = y1
           a(16,k) = 1.
        enddo

        call invert_matrix(n16, a, b)

        do k=1,n16
           van2(k,i,j) = b(k,n16)
        enddo

        sum0 = 0.
        do k=1,n16
           sum0 = sum0 + b(k,n16)
#ifdef CHECK_VAN2
           if ( k==1 .and. i==3 .and. j==3 ) then
                write(*,*) k,'Van2(3,3):', van2(k,i,j)
!               write(*,*) '          ', lon0, lat0
           endif
#endif
        enddo
        if (abs(sum0-1.)>1.e-12) then
            write(*,*) 'Failed van point:', i,j
            call mpp_error(FATAL, 'van2_init')
        endif
2000 continue
2500 continue


 end subroutine van2_init


#ifdef USE_EXTEND_CUBE
 subroutine extend_cube_s(non_ortho, grid, agrid, npx, npy, symm)
 
! Initialization of interpolation factors for the extended cubed sphere
! for interpolating cell mean scalars beyond the cubed face
 
 logical, intent(in):: non_ortho
 real,    intent(in)::  grid(isd:ied+1,jsd:jed+1,2)
 real,    intent(in):: agrid(isd:ied  ,jsd:jed  ,2)
 integer, intent(in):: npx, npy
 logical, intent(in):: symm  ! Not working; requires global grids

 real p1(3), p2(3), p3(3), p4(3), p5(3), pp(3)
 real q1(2), q2(2)
 real d1, d2, d3
 integer i, j
 integer im2, jm2
 logical local_in, local_out
 real, parameter:: esl = 1.E-5

 allocate ( ex_s(npx) )
 allocate ( ex_n(npx) )
 allocate ( ex_w(npy) )
 allocate ( ex_e(npy) )


  if ( .not. non_ortho ) then
     ex_s = 0.
     ex_n = 0.
     ex_w = 0.
     ex_e = 0.
  else
     ex_s = big_number 
     ex_n = big_number
     ex_w = big_number
     ex_e = big_number
 
     if ( npx /= npy ) call mpp_error(FATAL, 'extend_cube_s: npx /= npy')
     if ( (npx/2)*2 == npx ) call mpp_error(FATAL, 'extend_cube_s: npx/npy is not an odd number')

     im2 = (npx-1)/2
     jm2 = (npy-1)/2

 if ( is==1 ) then
    i=1
    do j=js,je
         call latlon2xyz( agrid(i,  j,  1:2), p1)
         call mid_pt_cart(grid(i,j,1:2), grid(i,j+1,1:2), p2)
       if ( j<=jm2 ) then
! q_w(j) = (1.-ex_w(j)) * q(j) + ex_w(j) * q(j+1)
! 1st column
          call latlon2xyz( agrid(i-1,j,  1:2), p3)
          call latlon2xyz( agrid(i-1,j+1,1:2), p4)
          call intersect(p1, p2, p3, p4, 1., pp, local_in, local_out)
          call cart_to_latlon(1, pp, q1(1), q1(2))
          d1 = great_circle_dist( q1, agrid(i-1,j,  1:2) )
          d2 = great_circle_dist( q1, agrid(i-1,j+1,1:2) )
          d3 = great_circle_dist( agrid(i-1,j,1:2), agrid(i-1,j+1,1:2) )
          if ( d1 > d3 ) then
               call mpp_error(FATAL, 'extend_cube_s: 1st column intp violated')
          else
               ex_w(j) = d1 / ( d1 + d2 )
          endif
          if( ex_w(j) < esl ) ex_w(j) = 0.
!         if(gid==0) write(*,*) i,j, ex_w(j)
       else
!
! q_w(j) = (1.-ex_w(j)) * q(j) + ex_w(j) * q(j-1)
! 1st column
          call latlon2xyz( agrid(i-1,j,  1:2), p3)
          call latlon2xyz( agrid(i-1,j-1,1:2), p4)
          call intersect(p1, p2, p3, p4, 1., pp, local_in, local_out)
          call cart_to_latlon(1, pp, q1(1), q1(2))
          d1 = great_circle_dist( q1, agrid(i-1,j,  1:2) )
          d2 = great_circle_dist( q1, agrid(i-1,j-1,1:2) )
          d3 = great_circle_dist( agrid(i-1,j,1:2), agrid(i-1,j-1,1:2) )
          if ( d1 > d3 ) then
               call mpp_error(FATAL, 'extend_cube_s: 1st column intp violated')
          else
               ex_w(j) = d1 / ( d1 + d2 )
          endif
          if( ex_w(j) < esl ) ex_w(j) = 0.
!         if(gid==0) write(*,*) i,j, ex_w(j)
       endif
    enddo
 endif

 if ( (ie+1)==npx ) then
    i=npx-1
    do j=js,je
         call latlon2xyz( agrid(i  ,j, 1:2), p1)
         call mid_pt_cart(grid(i+1,j,1:2), grid(i+1,j+1,1:2), p2)
       if ( j<=jm2 ) then
! q_e(j) = (1.-ex_e(j)) * q(j) + ex_e(j) * q(j+1)
! 1st column
          call latlon2xyz( agrid(i+1,j,  1:2), p3)
          call latlon2xyz( agrid(i+1,j+1,1:2), p4)
          call intersect(p1, p2, p3, p4, 1., pp, local_in, local_out)
          call cart_to_latlon(1, pp, q1(1), q1(2))
          d1 = great_circle_dist( q1, agrid(i+1,j,  1:2) )
          d2 = great_circle_dist( q1, agrid(i+1,j+1,1:2) )
          d3 = great_circle_dist( agrid(i+1,j,1:2), agrid(i+1,j+1,1:2) )
          if ( d1 > d3 ) then
               call mpp_error(FATAL, 'extend_cube_s: 1st column intp violated')
          else
               ex_e(j) = d1 / ( d1 + d2 )
          endif
          if( ex_e(j) < esl ) ex_e(j) = 0.
!         if(gid==0) write(*,*) i,j, ex_e(j) - ex_w(j)
       else
!
! q_e(j) = (1.-ex_e(j)) * q(j) + ex_e(j) * q(j-1)
! 1st column
          call latlon2xyz( agrid(i+1,j,  1:2), p3)
          call latlon2xyz( agrid(i+1,j-1,1:2), p4)
          call intersect(p1, p2, p3, p4, 1., pp, local_in, local_out)
          call cart_to_latlon(1, pp, q1(1), q1(2))
          d1 = great_circle_dist( q1, agrid(i+1,j,  1:2) )
          d2 = great_circle_dist( q1, agrid(i+1,j-1,1:2) )
          d3 = great_circle_dist( agrid(i+1,j,1:2), agrid(i+1,j-1,1:2) )
          if ( d1 > d3 ) then
               call mpp_error(FATAL, 'extend_cube_s: 1st column intp violated')
          else
               ex_e(j) = d1 / ( d1 + d2 )
          endif
          if( ex_e(j) < esl ) ex_e(j) = 0.
!         if(gid==0) write(*,*) i,j, ex_e(j) - ex_w(j)
       endif
    enddo
 endif

! Make it symmetrical
 if ( symm) then
    do j=js,je
       ex_e(j) = 0.5*(ex_e(j) + ex_w(j))
       ex_w(j) = ex_e(j)
    enddo
 endif

 if ( js==1 ) then
    j=1
    do i=is,ie
          call latlon2xyz( agrid(i,j,  1:2), p1)
          call mid_pt_cart(grid(i,j,1:2), grid(i+1,j,1:2), p2)
       if ( i<=im2 ) then
! q_s(i) = (1.-ex_s(i)) * q(i) + ex_s(i) * q(i+1)
! 1st row
          call latlon2xyz( agrid(i,  j-1,1:2), p3)
          call latlon2xyz( agrid(i+1,j-1,1:2), p4)
          call intersect(p1, p2, p3, p4, 1., pp, local_in, local_out)
          call cart_to_latlon(1, pp, q1(1), q1(2))
          d1 = great_circle_dist( q1, agrid(i,  j-1,1:2) )
          d2 = great_circle_dist( q1, agrid(i+1,j-1,1:2) )
          d3 = great_circle_dist( agrid(i,j-1,1:2), agrid(i+1,j-1,1:2) )
          if ( d1 > d3 ) then
               call mpp_error(FATAL, 'extend_cube_s: 1st column intp violated')
          else
               ex_s(i) = d1 / ( d1 + d2 )
          endif
          if( ex_s(i) < esl ) ex_s(i) = 0.
!         if(gid==0) write(*,*) i,j, ex_s(i)
       else
! q_s(i) = (1.-ex_s(i)) * q(i) + ex_s(i) * q(i-1)
! 1st row
          call latlon2xyz( agrid(i,  j-1,1:2), p3)
          call latlon2xyz( agrid(i-1,j-1,1:2), p4)
          call intersect(p1, p2, p3, p4, 1., pp, local_in, local_out)
          call cart_to_latlon(1, pp, q1(1), q1(2))
          d1 = great_circle_dist( q1, agrid(i,  j-1,1:2) )
          d2 = great_circle_dist( q1, agrid(i-1,j-1,1:2) )
          d3 = great_circle_dist( agrid(i,j-1,1:2), agrid(i-1,j-1,1:2) )
          if ( d1 > d3 ) then
               call mpp_error(FATAL, 'extend_cube_s: 1st column intp violated')
          else
               ex_s(i) = d1 / ( d1 + d2 )
          endif
          if( ex_s(i) < esl ) ex_s(i) = 0.
!         if(gid==0) write(*,*) i,j, ex_s(i)
       endif
    enddo
 endif


 if ( (je+1)==npy ) then
    j=npy-1
    do i=is,ie
          call latlon2xyz( agrid(i,j,  1:2), p1)
          call mid_pt_cart(grid(i,j+1,1:2), grid(i+1,j+1,1:2), p2)
       if ( i<=im2 ) then
! q_n(i) = (1.-ex_n(i)) * q(i) + ex_n(i) * q(i+1)
! 1st row
          call latlon2xyz( agrid(i,  j+1,1:2), p3)
          call latlon2xyz( agrid(i+1,j+1,1:2), p4)
          call intersect(p1, p2, p3, p4, 1., pp, local_in, local_out)
          call cart_to_latlon(1, pp, q1(1), q1(2))
          d1 = great_circle_dist( q1, agrid(i,  j+1,1:2) )
          d2 = great_circle_dist( q1, agrid(i+1,j+1,1:2) )
          d3 = great_circle_dist( agrid(i,j+1,1:2), agrid(i+1,j+1,1:2) )
          if ( d1 > d3 ) then
               call mpp_error(FATAL, 'extend_cube_s: 1st column intp violated')
          else
               ex_n(i) = d1 / ( d1 + d2 )
          endif
          if( ex_n(i) < esl ) ex_n(i) = 0.
!         if(gid==0) write(*,*) i,j, ex_n(i) - ex_s(i)
       else
! q_n(i) = (1.-ex_n(i)) * q(i) + ex_n(i) * q(i-1)
! 1st row
          call latlon2xyz( agrid(i,  j+1,1:2), p3)
          call latlon2xyz( agrid(i-1,j+1,1:2), p4)
          call intersect(p1, p2, p3, p4, 1., pp, local_in, local_out)
          call cart_to_latlon(1, pp, q1(1), q1(2))
          d1 = great_circle_dist( q1, agrid(i,  j+1,1:2) )
          d2 = great_circle_dist( q1, agrid(i-1,j+1,1:2) )
          d3 = great_circle_dist( agrid(i,j+1,1:2), agrid(i-1,j+1,1:2) )
          if ( d1 > d3 ) then
               call mpp_error(FATAL, 'extend_cube_s: 1st column intp violated')
          else
               ex_n(i) = d1 / ( d1 + d2 )
          endif
          if( ex_n(i) < esl ) ex_n(i) = 0.
!         if(gid==0) write(*,*) i,j, ex_n(i) - ex_s(i)
       endif
    enddo
 endif

! Make it symmetrical
 if ( symm) then
    do i=is,ie
       ex_n(i) = 0.5*(ex_n(i) + ex_s(i))
       ex_s(i) = ex_n(i)
    enddo
 endif

 endif

 end subroutine extend_cube_s
#endif


 subroutine efactor_a2c_v(non_ortho, grid, agrid, npx, npy)
!
! Initialization of interpolation factors at face edges
! for interpolating vectors from A to C grid
!
 logical, intent(in):: non_ortho
 real,    intent(in)::  grid(isd:ied+1,jsd:jed+1,2)
 real,    intent(in):: agrid(isd:ied  ,jsd:jed  ,2)
 integer, intent(in):: npx, npy

 real px(2,isd:ied+1),  py(2,jsd:jed+1)
 real p1(2,isd:ied+1),  p2(2,jsd:jed+1)       ! mid-point
 real d1, d2
 integer i, j
 integer im2, jm2

 allocate ( edge_vect_s(isd:ied) )
 allocate ( edge_vect_n(isd:ied) )
 allocate ( edge_vect_w(jsd:jed) )
 allocate ( edge_vect_e(jsd:jed) )

  if ( .not. non_ortho ) then
     edge_vect_s = 0.
     edge_vect_n = 0.
     edge_vect_w = 0.
     edge_vect_e = 0.
  else
     edge_vect_s = big_number
     edge_vect_n = big_number
     edge_vect_w = big_number
     edge_vect_e = big_number

     if ( npx /= npy ) call mpp_error(FATAL, 'efactor_a2c_v: npx /= npy')
     if ( (npx/2)*2 == npx ) call mpp_error(FATAL, 'efactor_a2c_v: npx/npy is not an odd number')

     im2 = (npx-1)/2
     jm2 = (npy-1)/2

 if ( is==1 ) then
    i=1
    do j=js-2,je+2
       call mid_pt_sphere(agrid(i-1,j,1:2), agrid(i,j,  1:2), py(1,j))
       call mid_pt_sphere( grid(i,  j,1:2),  grid(i,j+1,1:2), p2(1,j))
    enddo

! west edge:
!------------------------------------------------------------------
! v_sw(j) = (1.-edge_vect_w(j)) * p(j) + edge_vect_w(j) * p(j+1)
!------------------------------------------------------------------
    do j=js-1,je+1
       if ( j<=jm2 ) then
            d1 = great_circle_dist( py(1,j  ), p2(1,j) )
            d2 = great_circle_dist( py(1,j+1), p2(1,j) )
            edge_vect_w(j) = d1 / ( d1 + d2 )
       else
            d2 = great_circle_dist( py(1,j-1), p2(1,j) )
            d1 = great_circle_dist( py(1,j  ), p2(1,j) )
            edge_vect_w(j) = d1 / ( d2 + d1 )
       endif
    enddo
    if ( js==1 ) then
         edge_vect_w(0) = edge_vect_w(1)
    endif
    if ( (je+1)==npy ) then
         edge_vect_w(npy) = edge_vect_w(je)
    endif
    do j=js-1,je+1
!      if ( gid==0 ) write(*,*) j, edge_vect_w(j)
    enddo
 endif

 if ( (ie+1)==npx ) then
    i=npx
    do j=jsd,jed
       call mid_pt_sphere(agrid(i-1,j,1:2), agrid(i,j,  1:2), py(1,j))
       call mid_pt_sphere( grid(i,  j,1:2),  grid(i,j+1,1:2), p2(1,j))
    enddo

    do j=js-1,je+1
       if ( j<=jm2 ) then
            d1 = great_circle_dist( py(1,j  ), p2(1,j) )
            d2 = great_circle_dist( py(1,j+1), p2(1,j) )
            edge_vect_e(j) = d1 / ( d1 + d2 )
       else
            d2 = great_circle_dist( py(1,j-1), p2(1,j) )
            d1 = great_circle_dist( py(1,j  ), p2(1,j) )
            edge_vect_e(j) = d1 / ( d2 + d1 )
       endif
    enddo
    if ( js==1 ) then
         edge_vect_e(0) = edge_vect_e(1)
    endif
    if ( (je+1)==npy ) then
         edge_vect_e(npy) = edge_vect_e(je)
    endif
    do j=js-1,je+1
!      if ( gid==0 ) write(*,*) j, edge_vect_e(j)
    enddo
 endif

 if ( js==1 ) then
    j=1
    do i=isd,ied
       call mid_pt_sphere(agrid(i,j-1,1:2), agrid(i,  j,1:2), px(1,i))
       call mid_pt_sphere( grid(i,j,  1:2),  grid(i+1,j,1:2), p1(1,i))
    enddo
! south_west edge:
!------------------------------------------------------------------
! v_s(i) = (1.-edge_vect_s(i)) * p(i) + edge_vect_s(i) * p(i+1)
!------------------------------------------------------------------
    do i=is-1,ie+1
       if ( i<=im2 ) then
            d1 = great_circle_dist( px(1,i  ), p1(1,i) )
            d2 = great_circle_dist( px(1,i+1), p1(1,i) )
            edge_vect_s(i) = d1 / ( d1 + d2 )
       else
            d2 = great_circle_dist( px(1,i-1), p1(1,i) )
            d1 = great_circle_dist( px(1,i  ), p1(1,i) )
            edge_vect_s(i) = d1 / ( d2 + d1 )
       endif
    enddo
    if ( is==1 ) then
         edge_vect_s(0) = edge_vect_s(1)
    endif
    if ( (ie+1)==npx ) then
         edge_vect_s(npx) = edge_vect_s(ie)
    endif
    do i=is-1,ie+1
!      if ( gid==0 ) write(*,*) i, edge_vect_s(i)
    enddo
 endif


 if ( (je+1)==npy ) then
! v_n(i) = (1.-edge_vect_n(i)) * p(i) + edge_vect_n(i) * p(i+1)
    j=npy
    do i=isd,ied
       call mid_pt_sphere(agrid(i,j-1,1:2), agrid(i,  j,1:2), px(1,i))
       call mid_pt_sphere( grid(i,j,  1:2),  grid(i+1,j,1:2), p1(1,i))
    enddo

    do i=is-1,ie+1
       if ( i<=im2 ) then
            d1 = great_circle_dist( px(1,i  ), p1(1,i) )
            d2 = great_circle_dist( px(1,i+1), p1(1,i) )
            edge_vect_n(i) = d1 / ( d1 + d2 )
       else
            d2 = great_circle_dist( px(1,i-1), p1(1,i) )
            d1 = great_circle_dist( px(1,i  ), p1(1,i) )
            edge_vect_n(i) = d1 / ( d2 + d1 )
       endif
    enddo
    if ( is==1 ) then
         edge_vect_n(0) = edge_vect_n(1)
    endif
    if ( (ie+1)==npx ) then
         edge_vect_n(npx) = edge_vect_n(ie)
    endif
    do i=is-1,ie+1
!      if ( gid==0 ) write(*,*) i, edge_vect_n(i)
    enddo
 endif

 endif

 end subroutine efactor_a2c_v


 subroutine edge_factors(non_ortho, grid, agrid, npx, npy)
!
! Initialization of interpolation factors at face edges
! for interpolation from A to B grid
!
 logical, intent(in):: non_ortho
 real,    intent(in)::  grid(isd:ied+1,jsd:jed+1,2)
 real,    intent(in):: agrid(isd:ied  ,jsd:jed  ,2)
 integer, intent(in):: npx, npy

 real px(2,npx), py(2,npy)
 real d1, d2
 integer i, j

 allocate ( edge_s(npx) )
 allocate ( edge_n(npx) )
 allocate ( edge_w(npy) )
 allocate ( edge_e(npy) )


  if ( .not. non_ortho ) then
     edge_s = 0.5
     edge_n = 0.5
     edge_w = 0.5
     edge_e = 0.5
  else
     edge_s = big_number
     edge_n = big_number
     edge_w = big_number
     edge_e = big_number
 
! west edge:
!----------------------------------------------------------
! p_west(j) = (1.-edge_w(j)) * p(j) + edge_w(j) * p(j-1)
!----------------------------------------------------------
 if ( is==1 ) then
    i=1
    do j=max(1,js-1), min(npy-1,je+1)
       call mid_pt_sphere(agrid(i-1,j,1:2), agrid(i,j,1:2), py(1,j))
    enddo
    do j=max(2,js), min(npy-1,je+1)
       d1 = great_circle_dist( py(1,j-1), grid(i,j,1:2) )
       d2 = great_circle_dist( py(1,j  ), grid(i,j,1:2) )
       edge_w(j) = d2 / ( d1 + d2 )
    enddo
 endif

! east edge:
!----------------------------------------------------------
! p_east(j) = (1.-edge_e(j)) * p(j) + edge_e(j) * p(j-1)
!----------------------------------------------------------
 if ( (ie+1)==npx ) then
    i=npx
    do j=max(1,js-1), min(npy-1,je+1)
       call mid_pt_sphere(agrid(i-1,j,1:2), agrid(i,j,1:2), py(1,j))
    enddo
    do j=max(2,js), min(npy-1,je+1)
       d1 = great_circle_dist( py(1,j-1), grid(i,j,1:2) )
       d2 = great_circle_dist( py(1,j  ), grid(i,j,1:2) )
       edge_e(j) = d2 / ( d1 + d2 )
! Check rounding difference:
!      if(gid==0) write(*,*) j, edge_w(j) - edge_e(j)
    enddo
 endif


! south edge:
!----------------------------------------------------------
! p_south(j) = (1.-edge_s(i)) * p(i) + edge_s(i) * p(i-1)
!----------------------------------------------------------
 if ( js==1 ) then
    j=1
    do i=max(1,is-1), min(npx-1,ie+1)
       call mid_pt_sphere(agrid(i,j-1,1:2), agrid(i,j,1:2), px(1,i))
    enddo
    do i=max(2,is), min(npx-1,ie+1)
       d1 = great_circle_dist( px(1,i-1), grid(i,j,1:2) )
       d2 = great_circle_dist( px(1,i  ), grid(i,j,1:2) )
       edge_s(i) = d2 / ( d1 + d2 )
    enddo
 endif

! North edge:
!----------------------------------------------------------
! p_north(j) = (1.-edge_n(i)) * p(i) + edge_n(i) * p(i-1)
!----------------------------------------------------------
 if ( (je+1)==npy ) then
    j=npy
    do i=max(1,is-1), min(npx-1,ie+1)
       call mid_pt_sphere(agrid(i,j-1,1:2), agrid(i,j,1:2), px(1,i))
    enddo
    do i=max(2,is), min(npx-1,ie+1)
       d1 = great_circle_dist( px(1,i-1), grid(i,j,1:2) )
       d2 = great_circle_dist( px(1,i  ), grid(i,j,1:2) )
       edge_n(i) = d2 / ( d1 + d2 )
!      if(gid==0) write(*,*) i, edge_s(i), edge_n(i)-edge_s(i)
    enddo
 endif
 endif

 end subroutine edge_factors


 subroutine gnomonic_grids(grid_type, im, lon, lat)
 integer, intent(in):: im, grid_type
 real, intent(out):: lon(im+1,im+1)
 real, intent(out):: lat(im+1,im+1)
 real pi
 integer i, j

  pi = 4.*atan(1.)


  if(grid_type==0) call gnomonic_ed(  im, lon, lat)
  if(grid_type==1) call gnomonic_dist(im, lon, lat)
  if(grid_type==2) call gnomonic_angl(im, lon, lat)


  if(grid_type<3) then
     call symm_ed(im, lon, lat)
     do j=1,im+1
        do i=1,im+1
           lon(i,j) = lon(i,j) - pi
        enddo
     enddo
!    call van2_init(lon, lat, im+1, im+1)
  endif
  
 end subroutine gnomonic_grids



 subroutine gnomonic_ed(im, lamda, theta)
!-----------------------------------------------------
! Equal distance along the 4 edges of the cubed sphere
!-----------------------------------------------------
! Properties: 
!            * defined by intersections of great circles
!            * max(dx,dy; global) / min(dx,dy; global) = sqrt(2) = 1.4142
!            * Max(aspect ratio) = 1.06089
!            * the N-S coordinate curves are const longitude on the 4 faces with equator 
! For C2000: (dx_min, dx_max) = (3.921, 5.545)    in km unit
! This is the grid of choice for global cloud resolving

 integer, intent(in):: im
 real, intent(out):: lamda(im+1,im+1)
 real, intent(out):: theta(im+1,im+1)

! Local:
 real pp(3,im+1,im+1)
! real(f_p):: pi, rsq3, alpha, delx, dely
 real:: pi, rsq3, alpha, delx, dely
 integer i, j, k

    pi = 4.*atan(1.)
  rsq3 = 1./sqrt(3.) 
 alpha = asin( rsq3 )

! Ranges:
! lamda = [0.75*pi, 1.25*pi]
! theta = [-alpha, alpha]

    dely = 2.*alpha / real(im)

! Define East-West edges:
 do j=1,im+1
    lamda(1,   j) = 0.75*pi                  ! West edge
    lamda(im+1,j) = 1.25*pi                  ! East edge
    theta(1,   j) = -alpha + dely*real(j-1)  ! West edge
    theta(im+1,j) = theta(1,j)               ! East edge
 enddo

! Get North-South edges by symmetry:

 do i=2,im
    call mirror_latlon(lamda(1,1), theta(1,1), lamda(im+1,im+1), theta(im+1,im+1), &
                       lamda(1,i), theta(1,i), lamda(i,1),       theta(i,      1) )
    lamda(i,im+1) =  lamda(i,1)
    theta(i,im+1) = -theta(i,1)
 enddo

! Set 4 corners:
    call latlon2xyz2(lamda(1    ,  1), theta(1,      1), pp(1,   1,   1))
    call latlon2xyz2(lamda(im+1,   1), theta(im+1,   1), pp(1,im+1,   1))
    call latlon2xyz2(lamda(1,   im+1), theta(1,   im+1), pp(1,   1,im+1))
    call latlon2xyz2(lamda(im+1,im+1), theta(im+1,im+1), pp(1,im+1,im+1))

! Map edges on the sphere back to cube:
! Intersections at x=-rsq3

 i=1
 do j=2,im
    call latlon2xyz2(lamda(i,j), theta(i,j), pp(1,i,j))
    pp(2,i,j) = -pp(2,i,j)*rsq3/pp(1,i,j)
    pp(3,i,j) = -pp(3,i,j)*rsq3/pp(1,i,j)
 enddo

 j=1
 do i=2,im
    call latlon2xyz2(lamda(i,j), theta(i,j), pp(1,i,1))
    pp(2,i,1) = -pp(2,i,1)*rsq3/pp(1,i,1)
    pp(3,i,1) = -pp(3,i,1)*rsq3/pp(1,i,1)
 enddo

 do j=1,im+1
    do i=1,im+1
       pp(1,i,j) = -rsq3
    enddo
 enddo

 do j=2,im+1
    do i=2,im+1
! Copy y-z face of the cube along j=1
       pp(2,i,j) = pp(2,i,1)
! Copy along i=1
       pp(3,i,j) = pp(3,1,j)
    enddo
 enddo

 call cart_to_latlon( (im+1)*(im+1), pp, lamda, theta)

 end subroutine gnomonic_ed



 subroutine gnomonic_angl(im, lamda, theta)
! This is the commonly known equi-angular grid
 integer im
 real lamda(im+1,im+1)
 real theta(im+1,im+1)
 real p(3,im+1,im+1)
! Local
 real rsq3, xf, y0, z0, y, x, z, ds
 real dy, dz
 integer j,k
 real pi, dp

 pi = 4.*atan(1.)
 dp = 0.5*pi/real(im)

 rsq3 = 1./sqrt(3.) 
 do k=1,im+1
    do j=1,im+1
       p(1,j,k) =-rsq3               ! constant
       p(2,j,k) =-rsq3*tan(-0.25*pi+(j-1)*dp)
       p(3,j,k) = rsq3*tan(-0.25*pi+(k-1)*dp)
    enddo
 enddo

 call cart_to_latlon( (im+1)*(im+1), p, lamda, theta)

 end subroutine gnomonic_angl

 subroutine gnomonic_dist(im, lamda, theta)
! This is the commonly known equi-distance grid
 integer im
 real lamda(im+1,im+1)
 real theta(im+1,im+1)
 real p(3,im+1,im+1)
! Local
 real rsq3, xf, y0, z0, y, x, z, ds
 real dy, dz
 integer j,k
 real pi

 pi = 4.*atan(1.)

! Face-2

 rsq3 = 1./sqrt(3.) 
 xf = -rsq3
 y0 =  rsq3;  dy = -2.*rsq3/im 
 z0 = -rsq3;  dz =  2.*rsq3/im

 do k=1,im+1
    do j=1,im+1
       p(1,j,k) = xf
       p(2,j,k) = y0 + (j-1)*dy
       p(3,j,k) = z0 + (k-1)*dz
    enddo
 enddo
 call cart_to_latlon( (im+1)*(im+1), p, lamda, theta)

 end subroutine gnomonic_dist

 subroutine symm_ed(im, lamda, theta)
! Make grid symmetrical to i=im/2+1
 integer im
 real lamda(im+1,im+1)
 real theta(im+1,im+1)
 integer i,j,ip,jp
 real pi, avg

 pi = 4.*atan(1.)

 do j=2,im+1
    do i=2,im
       lamda(i,j) = lamda(i,1)
    enddo
 enddo

 do j=1,im+1
    do i=1,im/2
       ip = im + 2 - i
       avg = 0.5*(lamda(i,j)-lamda(ip,j))
       lamda(i, j) = avg + pi
       lamda(ip,j) = pi - avg 
       avg = 0.5*(theta(i,j)+theta(ip,j))
       theta(i, j) = avg
       theta(ip,j) = avg
    enddo
 enddo

! Make grid symmetrical to j=im/2+1
 do j=1,im/2
       jp = im + 2 - j
    do i=2,im
       avg = 0.5*(lamda(i,j)+lamda(i,jp))
       lamda(i, j) = avg
       lamda(i,jp) = avg
       avg = 0.5*(theta(i,j)-theta(i,jp))
       theta(i, j) =  avg
       theta(i,jp) = -avg
    enddo
 enddo

 end subroutine symm_ed

 subroutine latlon2xyz2(lon, lat, p3)
 real, intent(in):: lon, lat
 real, intent(out):: p3(3)
 real e(2)

    e(1) = lon;    e(2) = lat
    call latlon2xyz(e, p3)

 end subroutine latlon2xyz2


 subroutine latlon2xyz(p, e)
!
! Routine to map (lon, lat) to (x,y,z)
!
 real, intent(in) :: p(2)
 real, intent(out):: e(3)

 integer n
 real (f_p):: q(2)
 real (f_p):: e1, e2, e3

    do n=1,2
       q(n) = p(n)
    enddo

    e1 = cos(q(2)) * cos(q(1))
    e2 = cos(q(2)) * sin(q(1))
    e3 = sin(q(2))
!-----------------------------------
! Truncate to the desired precision:
!-----------------------------------
    e(1) = e1
    e(2) = e2
    e(3) = e3

 end subroutine latlon2xyz


 subroutine mirror_xyz(p1, p2, p0, p)

! Given the "mirror" as defined by p1(x1, y1, z1), p2(x2, y2, z2), and center 
! of the sphere, compute the mirror image of p0(x0, y0, z0) as p(x, y, z)

 real, intent(in) :: p1(3), p2(3), p0(3)
 real, intent(out):: p(3)
!
 real:: x1, y1, z1, x2, y2, z2, x0, y0, z0
 real nb(3)
 real pdot
 integer k

 call vect_cross(nb, p1, p2)
    pdot = sqrt(nb(1)**2+nb(2)**2+nb(3)**2)
 do k=1,3
    nb(k) = nb(k) / pdot
 enddo

 pdot = p0(1)*nb(1) + p0(2)*nb(2) + p0(3)*nb(3)
 do k=1,3
    p(k) = p0(k) - 2.*pdot*nb(k)
 enddo

 end subroutine mirror_xyz 


 subroutine mirror_latlon(lon1, lat1, lon2, lat2, lon0, lat0, lon3, lat3)
!
! Given the "mirror" as defined by (lon1, lat1), (lon2, lat2), and center 
! of the sphere, compute the mirror image of (lon0, lat0) as  (lon3, lat3)

 real, intent(in):: lon1, lat1, lon2, lat2, lon0, lat0
 real, intent(out):: lon3, lat3
!
 real p0(3), p1(3), p2(3), nb(3), pp(3), sp(2)
 real pdot
 integer k

 call latlon2xyz2(lon0, lat0, p0)
 call latlon2xyz2(lon1, lat1, p1)
 call latlon2xyz2(lon2, lat2, p2)
 call vect_cross(nb, p1, p2)

 pdot = sqrt(nb(1)**2+nb(2)**2+nb(3)**2)
 do k=1,3
    nb(k) = nb(k) / pdot
 enddo

 pdot = p0(1)*nb(1) + p0(2)*nb(2) + p0(3)*nb(3)
 do k=1,3
    pp(k) = p0(k) - 2.*pdot*nb(k)
 enddo

 call cart_to_latlon(1, pp, sp(1), sp(2))
 lon3 = sp(1)
 lat3 = sp(2)

 end subroutine  mirror_latlon


 subroutine cart_to_latlon(np, q, xs, ys)
! vector version of cart_to_latlon1
  integer, intent(in):: np
  real, intent(inout):: q(3,np)
  real, intent(inout):: xs(np), ys(np)
! local
  real, parameter:: esl=1.e-10
  real (f_p):: p(3)
  real (f_p):: pi, dist, lat, lon
  integer i,k

  pi = 4.*atan(1.)

  do i=1,np
     do k=1,3
        p(k) = q(k,i)
     enddo
     dist = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
     do k=1,3
        p(k) = p(k) / dist
     enddo

     if ( (abs(p(1))+abs(p(2)))  < esl ) then
          lon = 0.
     else
          lon = atan2( p(2), p(1) )   ! range [-pi,pi]
     endif

     if ( lon < 0.) lon = 2.*pi + lon
     lat = asin(p(3))
     
     xs(i) = lon
     ys(i) = lat
! q Normalized:
     do k=1,3
        q(k,i) = p(k)
     enddo
  enddo

 end  subroutine cart_to_latlon



 subroutine vect_cross(e, p1, p2)
 real, intent(in) :: p1(3), p2(3)
 real, intent(out):: e(3)
!
! Perform cross products of 3D vectors: e = P1 X P2
!
      e(1) = p1(2)*p2(3) - p1(3)*p2(2)
      e(2) = p1(3)*p2(1) - p1(1)*p2(3)
      e(3) = p1(1)*p2(2) - p1(2)*p2(1)

 end subroutine vect_cross



 subroutine get_center_vect( npx, npy, pp, u1, u2 )
    integer, intent(in):: npx, npy
    real, intent(in) :: pp(3,isd:ied+1,jsd:jed+1)
    real, intent(out):: u1(3,isd:ied,  jsd:jed)
    real, intent(out):: u2(3,isd:ied,  jsd:jed)
! Local:
    integer i,j,k
    real p1(3), p2(3), pc(3), p3(3)

    do j=jsd,jed
       do i=isd,ied
        if ( (i<1       .and. j<1  )     .or. (i>(npx-1) .and. j<1) .or.  &
             (i>(npx-1) .and. j>(npy-1)) .or. (i<1       .and. j>(npy-1))) then
             u1(1:3,i,j) = 0.
             u2(1:3,i,j) = 0.
        else
#ifdef NEW_VECT
          call cell_center3(pp(1,i,j), pp(1,i+1,j), pp(1,i,j+1), pp(1,i+1,j+1), pc)
! e1:
          call mid_pt3_cart(pp(1,i,j),   pp(1,i,j+1),   p1)
          call mid_pt3_cart(pp(1,i+1,j), pp(1,i+1,j+1), p2)
          call vect_cross(p3, p2, p1)
          call vect_cross(u1(1,i,j), pc, p3)
          call normalize_vect( u1(1,i,j) )
! e2:
          call mid_pt3_cart(pp(1,i,j),   pp(1,i+1,j),   p1)
          call mid_pt3_cart(pp(1,i,j+1), pp(1,i+1,j+1), p2)
          call vect_cross(p3, p2, p1)
          call vect_cross(u2(1,i,j), pc, p3)
          call normalize_vect( u2(1,i,j) )
#else
          do k=1,3
             u1(k,i,j) = pp(k,i+1,j)+pp(k,i+1,j+1) - pp(k,i,j)-pp(k,i,j+1)
             u2(k,i,j) = pp(k,i,j+1)+pp(k,i+1,j+1) - pp(k,i,j)-pp(k,i+1,j) 
          enddo
          call normalize_vect( u1(1,i,j) )
          call normalize_vect( u2(1,i,j) )
#endif
        endif
       enddo
    enddo

 end subroutine get_center_vect


 subroutine normalize_vect(e)
!                              Make e an unit vector
 real, intent(inout):: e(3)
 real(f_p):: pdot
 integer k

    pdot = e(1)**2 + e(2)**2 + e(3)**2
    pdot = sqrt( pdot ) 

!if ( pdot > 0. ) then
    do k=1,3
       e(k) = e(k) / pdot
    enddo
!else
!   do k=1,3
!      e(k) = 1. / sqrt(3.)
!   enddo
!endif

 end subroutine normalize_vect


 subroutine project_sphere_v( np, f, e )
!---------------------------------
 integer, intent(in):: np           ! total number of points
 real,    intent(in):: e(3,np)      ! input position unit vector
 real, intent(inout):: f(3,np)
! local
 real(f_p):: ap
 integer i

 do i=1,np
    ap = f(1,i)*e(1,i) + f(2,i)*e(2,i) + f(3,i)*e(3,i)
    f(1,i) = f(1,i) - ap*e(1,i)
    f(2,i) = f(2,i) - ap*e(2,i)
    f(3,i) = f(3,i) - ap*e(3,i)
 enddo

 end subroutine project_sphere_v


 subroutine intp_great_circle(beta, p1, p2, x_o, y_o)
 real, intent(in)::  beta    ! [0,1]
 real, intent(in)::  p1(2), p2(2)
 real, intent(out):: x_o, y_o     ! between p1 and p2 along GC
!------------------------------------------
    real:: pm(2)
    real:: e1(3), e2(3), e3(3)
    real:: s1, s2, s3, dd, alpha

      call latlon2xyz(p1, e1)
      call latlon2xyz(p2, e2)

       alpha = 1. - beta

       s1 = alpha*e1(1) + beta*e2(1)
       s2 = alpha*e1(2) + beta*e2(2)
       s3 = alpha*e1(3) + beta*e2(3)

       dd = sqrt( s1**2 + s2**2 + s3**2 )

       e3(1) = s1 / dd
       e3(2) = s2 / dd
       e3(3) = s3 / dd

      call cart_to_latlon(1, e3, pm(1), pm(2))

      x_o = pm(1)
      y_o = pm(2)

 end subroutine intp_great_circle


 subroutine mid_pt_sphere(p1, p2, pm)
      real , intent(IN)  :: p1(2), p2(2)
      real , intent(OUT) :: pm(2)
!------------------------------------------
      real e1(3), e2(3), e3(3)

      call latlon2xyz(p1, e1)
      call latlon2xyz(p2, e2)
      call mid_pt3_cart(e1, e2, e3)
      call cart_to_latlon(1, e3, pm(1), pm(2))

 end subroutine mid_pt_sphere



 subroutine mid_pt3_cart(p1, p2, e)
       real, intent(IN)  :: p1(3), p2(3)
       real, intent(OUT) :: e(3)
!
       real (f_p):: q1(3), q2(3)
       real (f_p):: dd, e1, e2, e3
       integer k

       do k=1,3
          q1(k) = p1(k)
          q2(k) = p2(k)
       enddo

       e1 = q1(1) + q2(1)
       e2 = q1(2) + q2(2)
       e3 = q1(3) + q2(3)

       dd = sqrt( e1**2 + e2**2 + e3**2 )
       e1 = e1 / dd
       e2 = e2 / dd
       e3 = e3 / dd

       e(1) = e1
       e(2) = e2
       e(3) = e3

 end subroutine mid_pt3_cart



 subroutine mid_pt_cart(p1, p2, e3)
    real, intent(IN)  :: p1(2), p2(2)
    real, intent(OUT) :: e3(3)
!-------------------------------------
    real e1(3), e2(3)

    call latlon2xyz(p1, e1)
    call latlon2xyz(p2, e2)
    call mid_pt3_cart(e1, e2, e3)

 end subroutine mid_pt_cart



 real function great_circle_dist( q1, q2, radius )
      real, intent(IN)           :: q1(2), q2(2)
      real, intent(IN), optional :: radius
 
      real (f_p):: p1(2), p2(2)
      real (f_p):: beta
      integer n

      do n=1,2
         p1(n) = q1(n)
         p2(n) = q2(n)
      enddo

      beta = asin( sqrt( sin((p1(2)-p2(2))/2.)**2 + cos(p1(2))*cos(p2(2))*   &
                         sin((p1(1)-p2(1))/2.)**2 ) ) * 2.

      if ( present(radius) ) then
           great_circle_dist = radius * beta
      else
           great_circle_dist = beta   ! Returns the angle
      endif

  end function great_circle_dist



 subroutine intersect(a1,a2,b1,b2,radius,x_inter,local_a,local_b)
  !--------------------------------------------------------------------!
  ! author:  Michael Herzog                                            !
  ! email:   Michael.Herzog@noaa.gov                                   !
  ! date:    July 2006                                                 !
  ! version: 0.1                                                       !
  !                                                                    !
  ! calculate intersection of two great circles                        !
  !--------------------------------------------------------------------!
    !------------------------------------------------------------------!
    ! calculate intersection of two great circles                      !
    !                                                                  !
    ! input:                                                           !
    ! a1, a2,  -   pairs of points on sphere in cartesian coordinates  !
    ! b1, b2       defining great circles                              !
    ! radius   -   radius of the sphere                                !
    !                                                                  !
    ! output:                                                          !
    ! x_inter  -   nearest intersection point of the great circles     !
    ! local_a  -   true if x1 between (a1, a2)                         !
    ! local_b  -   true if x1 between (b1, b2)                         !
    !------------------------------------------------------------------!
    real, dimension(3), intent(in)  :: a1, a2, b1, b2
    real, intent(in) :: radius
    real, dimension(3), intent(out) :: x_inter
    logical, intent(out) :: local_a,local_b
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real :: a2_xy, b1_xy, b2_xy, a2_xz, b1_xz, b2_xz,                   &
            b1_xyz, b2_xyz, length
    !------------------------------------------------------------------!
    ! calculate intersection point                                     !
    !------------------------------------------------------------------!
    a2_xy=a2(1)*a1(2)-a2(2)*a1(1)
    b1_xy=b1(1)*a1(2)-b1(2)*a1(1)
    b2_xy=b2(1)*a1(2)-b2(2)*a1(1)

    a2_xz=a2(1)*a1(3)-a2(3)*a1(1)
    b1_xz=b1(1)*a1(3)-b1(3)*a1(1)
    b2_xz=b2(1)*a1(3)-b2(3)*a1(1)

    b1_xyz=b1_xy*a2_xz-b1_xz*a2_xy
    b2_xyz=b2_xy*a2_xz-b2_xz*a2_xy

    if (b1_xyz==0.0) then
       x_inter(:)=b1(:)
    elseif (b2_xyz==0.0) then
       x_inter(:)=b2(:)
    else
       x_inter(:)=b2(:)-b1(:)*b2_xyz/b1_xyz
       length=sqrt(x_inter(1)*x_inter(1)+x_inter(2)*x_inter(2)+x_inter(3)*x_inter(3))
       x_inter(:)=radius/length*x_inter(:)
    endif
    !------------------------------------------------------------------!
    ! check if intersection is between pairs of points on sphere       !
    !------------------------------------------------------------------!
    call get_nearest()
    call check_local(a1,a2,local_a)
    call check_local(b1,b2,local_b)

  contains
    !------------------------------------------------------------------!
    subroutine get_nearest()
      real, dimension(3) :: center, dx
      real :: dist1,dist2

      center(:)=0.25*(a1(:)+a2(:)+b1(:)+b2(:))
      dx(:)=+x_inter(:)-center(:)
      dist1=dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3)
      dx(:)=-x_inter(:)-center(:)
      dist2=dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3)

      if (dist2<dist1) x_inter(:)=-x_inter(:)

    end subroutine get_nearest
    !------------------------------------------------------------------!
    subroutine check_local(x1,x2,local)
      real, dimension(3), intent(in) :: x1,x2
      logical, intent(out) :: local

      real, dimension(3) :: dx
      real :: dist, dist1, dist2

      dx(:)=x1(:)-x2(:)
      dist=dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3)
    
      dx(:)=x1(:)-x_inter(:)
      dist1=dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3)
      dx(:)=x2(:)-x_inter(:)
      dist2=dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3)

      if (dist1<=dist .and. dist2<=dist) then
         local=.true.
      else
         local=.false.
      endif
      
    end subroutine check_local
    !------------------------------------------------------------------!
  end subroutine intersect



  subroutine unit_vect_latlon(pp, elon, elat)
      real, intent(IN)  :: pp(2)
      real, intent(OUT) :: elon(3), elat(3)

      real (f_p):: lon, lat
      real (f_p):: sin_lon, cos_lon, sin_lat, cos_lat

      lon = pp(1)
      lat = pp(2)

      sin_lon = sin(lon)
      cos_lon = cos(lon)
      sin_lat = sin(lat)
      cos_lat = cos(lat)

      elon(1) = -sin_lon
      elon(2) =  cos_lon
      elon(3) =  0.

      elat(1) = -sin_lat*cos_lon
      elat(2) = -sin_lat*sin_lon
      elat(3) =  cos_lat

  end subroutine unit_vect_latlon



  real function v_prod(v1, v2)
  real v1(3), v2(3)

       v_prod = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

  end function v_prod



  subroutine init_cubed_to_latlon( agrid, grid_type, ord )

  real,    intent(in) :: agrid(isd:ied,jsd:jed,2)
  integer, intent(in) :: grid_type
  integer, intent(in) :: ord
  integer i, j

   g_type = grid_type
  c2l_ord = ord

  if ( g_type < 4 ) then

     allocate (  z11(is-1:ie+1,js-1:je+1) )
     allocate (  z12(is-1:ie+1,js-1:je+1) )
     allocate (  z21(is-1:ie+1,js-1:je+1) )
     allocate (  z22(is-1:ie+1,js-1:je+1) )

     allocate (  a11(is-1:ie+1,js-1:je+1) )
     allocate (  a12(is-1:ie+1,js-1:je+1) )
     allocate (  a21(is-1:ie+1,js-1:je+1) )
     allocate (  a22(is-1:ie+1,js-1:je+1) )
!     allocate ( vlon(is-1:ie+1,js-1:je+1,3) )
!     allocate ( vlat(is-1:ie+1,js-1:je+1,3) )
     allocate ( vlon(is-2:ie+2,js-2:je+2,3) )
     allocate ( vlat(is-2:ie+2,js-2:je+2,3) )

!     do j=js-1,je+1
!        do i=is-1,ie+1
     do j=js-2,je+2
        do i=is-2,ie+2
           call unit_vect_latlon(agrid(i,j,1:2), vlon(i,j,1:3), vlat(i,j,1:3))
        enddo
     enddo

     do j=js-1,je+1
        do i=is-1,ie+1
           z11(i,j) =  v_prod(ec1(1,i,j), vlon(i,j,1:3))
           z12(i,j) =  v_prod(ec1(1,i,j), vlat(i,j,1:3))
           z21(i,j) =  v_prod(ec2(1,i,j), vlon(i,j,1:3))
           z22(i,j) =  v_prod(ec2(1,i,j), vlat(i,j,1:3))
!-------------------------------------------------------------------------
           a11(i,j) =  0.5*v_prod(ec2(1,i,j), vlat(i,j,1:3)) / sina_s(i,j)
           a12(i,j) = -0.5*v_prod(ec1(1,i,j), vlat(i,j,1:3)) / sina_s(i,j)
           a21(i,j) = -0.5*v_prod(ec2(1,i,j), vlon(i,j,1:3)) / sina_s(i,j)
           a22(i,j) =  0.5*v_prod(ec1(1,i,j), vlon(i,j,1:3)) / sina_s(i,j)
        enddo
     enddo
  endif

  end subroutine init_cubed_to_latlon


 subroutine cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km, mode)
 integer, intent(in) :: km
 integer, intent(in), optional:: mode   ! update if present
 real, intent(in) :: dx(isd:ied,jsd:jed+1)
 real, intent(in) :: dy(isd:ied+1,jsd:jed)
 real, intent(in) ::rdxa(isd:ied,  jsd:jed)
 real, intent(in) ::rdya(isd:ied,  jsd:jed)
 real, intent(inout):: u(isd:ied,jsd:jed+1,km)
 real, intent(inout):: v(isd:ied+1,jsd:jed,km)
 real, intent(out):: ua(isd:ied, jsd:jed,km)
 real, intent(out):: va(isd:ied, jsd:jed,km)

 if ( c2l_ord == 2 ) then
      call c2l_ord2(u, v, ua, va, dx, dy, rdxa, rdya, km)
 else
      call c2l_ord4(u, v, ua, va, dx, dy, rdxa, rdya, km, mode)
 endif

 end subroutine cubed_to_latlon


 subroutine c2l_ord4(u, v, ua, va, dx, dy, rdxa, rdya, km, mode)

  integer, intent(in) :: km
  integer, intent(in), optional:: mode   ! update if present
  real, intent(in) ::  dx(isd:ied,jsd:jed+1)
  real, intent(in) ::  dy(isd:ied+1,jsd:jed)
  real, intent(in) ::rdxa(isd:ied,  jsd:jed)
  real, intent(in) ::rdya(isd:ied,  jsd:jed)
  real, intent(inout):: u(isd:ied,jsd:jed+1,km)
  real, intent(inout):: v(isd:ied+1,jsd:jed,km)
  real, intent(out)::  ua(isd:ied, jsd:jed,km)
  real, intent(out)::  va(isd:ied, jsd:jed,km)
! Local 
! 4-pt Lagrange interpolation
  real, parameter:: a1 =  0.5625
  real, parameter:: a2 = -0.0625
  real, parameter:: c1 =  1.125
  real, parameter:: c2 = -0.125
  real utmp(is:ie,  js:je+1)
  real vtmp(is:ie+1,js:je)
  real wu(is:ie,  js:je+1)
  real wv(is:ie+1,js:je)
  integer i, j, k

  if ( present(mode) ) then
                                   call timing_on('COMM_TOTAL')
       call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
                                  call timing_off('COMM_TOTAL')
  endif

 do k=1,km
   if ( g_type < 4 ) then
     do j=max(2,js),min(npyy-2,je)
        do i=max(2,is),min(npxx-2,ie)
           utmp(i,j) = c2*(u(i,j-1,k)+u(i,j+2,k)) + c1*(u(i,j,k)+u(i,j+1,k))
           vtmp(i,j) = c2*(v(i-1,j,k)+v(i+2,j,k)) + c1*(v(i,j,k)+v(i+1,j,k))
        enddo
     enddo

    if ( js==1 ) then
         do i=is,ie+1
            wv(i,1) = v(i,1,k)*dy(i,1)
         enddo
         do i=is,ie
            vtmp(i,1) = (wv(i,1) + wv(i+1,1)) * rdya(i,1)
            utmp(i,1) = (u(i,1,k)*dx(i,1) + u(i,2,k)*dx(i,2)) * rdxa(i,1)
         enddo
    endif

    if ( (je+1)==npyy ) then
         j = npyy-1
         do i=is,ie+1
            wv(i,j) = v(i,j,k)*dy(i,j)
         enddo
         do i=is,ie
            vtmp(i,j) = (wv(i,j) + wv(i+1,j)) * rdya(i,j)
            utmp(i,j) = (u(i,j,k)*dx(i,j) + u(i,j+1,k)*dx(i,j+1)) * rdxa(i,j)
         enddo
    endif

    if ( is==1 ) then
      i = 1
      do j=js,je
         wv(1,j) = v(1,j,k)*dy(1,j)
         wv(2,j) = v(2,j,k)*dy(2,j)
      enddo
      do j=js,je+1
         wu(i,j) = u(i,j,k)*dx(i,j)
      enddo
      do j=js,je
         utmp(i,j) = (wu(i,j) + wu(i,  j+1)) * rdxa(i,j)
         vtmp(i,j) = (wv(i,j) + wv(i+1,j  )) * rdya(i,j)
      enddo
    endif

    if ( (ie+1)==npxx ) then
      i = npxx-1
      do j=js,je
         wv(i,  j) = v(i,  j,k)*dy(i,  j)
         wv(i+1,j) = v(i+1,j,k)*dy(i+1,j)
      enddo
      do j=js,je+1
         wu(i,j) = u(i,j,k)*dx(i,j)
      enddo
      do j=js,je
         utmp(i,j) = (wu(i,j) + wu(i,  j+1)) * rdxa(i,j)
         vtmp(i,j) = (wv(i,j) + wv(i+1,j  )) * rdya(i,j)
      enddo
    endif

     do j=js,je
        do i=is,ie
           ua(i,j,k) = a11(i,j)*utmp(i,j) + a12(i,j)*vtmp(i,j)
           va(i,j,k) = a21(i,j)*utmp(i,j) + a22(i,j)*vtmp(i,j)
        enddo
     enddo
   else
! Simple Cartesian Geometry:
     do j=js,je
        do i=is,ie
           ua(i,j,k) = a2*(u(i,j-1,k)+u(i,j+2,k)) + a1*(u(i,j,k)+u(i,j+1,k))
           va(i,j,k) = a2*(v(i-1,j,k)+v(i+2,j,k)) + a1*(v(i,j,k)+v(i+1,j,k))
        enddo
     enddo
   endif
 enddo
 end subroutine c2l_ord4

 subroutine c2l_ord2(u, v, ua, va, dx, dy, rdxa, rdya, km)
  integer, intent(in) :: km
  real, intent(in) ::  u(isd:ied,jsd:jed+1,km)
  real, intent(in) ::  v(isd:ied+1,jsd:jed,km)
  real, intent(in) :: dx(isd:ied,jsd:jed+1)
  real, intent(in) :: dy(isd:ied+1,jsd:jed)
  real, intent(in) ::rdxa(isd:ied,  jsd:jed)
  real, intent(in) ::rdya(isd:ied,  jsd:jed)
!
  real, intent(out):: ua(isd:ied, jsd:jed,km)
  real, intent(out):: va(isd:ied, jsd:jed,km)
!--------------------------------------------------------------
! Local 
  real wu(is:ie,  js:je+1)
  real wv(is:ie+1,js:je)
  real u1(is:ie), v1(is:ie)
  integer i, j, k

!$omp parallel do default(shared) private(i, j, k, wu, wv, u1, v1)
  do k=1,km
     if ( g_type < 4 ) then
       do j=js,je+1
          do i=is,ie
             wu(i,j) = u(i,j,k)*dx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             wv(i,j) = v(i,j,k)*dy(i,j)
          enddo
       enddo

       do j=js,je
          do i=is,ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
             u1(i) = (wu(i,j) + wu(i,j+1)) * rdxa(i,j)
             v1(i) = (wv(i,j) + wv(i+1,j)) * rdya(i,j)
! Cubed (cell center co-variant winds) to lat-lon:
             ua(i,j,k) = a11(i,j)*u1(i) + a12(i,j)*v1(i)
             va(i,j,k) = a21(i,j)*u1(i) + a22(i,j)*v1(i)
          enddo
       enddo
     else
! 2nd order:
       do j=js,je
          do i=is,ie
             ua(i,j,k) = 0.5*(u(i,j,k)+u(i,  j+1,k))
             va(i,j,k) = 0.5*(v(i,j,k)+v(i+1,j,  k))
          enddo
       enddo
     endif
  enddo

 end subroutine c2l_ord2


 subroutine expand_cell(q1, q2, q3, q4, a1, a2, a3, a4, fac)
! Util for land model (for BW)
!
!        4----3
!        |  . |
!        1----2
!
      real, intent(in):: q1(2), q2(2), q3(2), q4(2)
      real, intent(in):: fac    ! expansion factor: outside: > 1
                                ! fac = 1: qq1 returns q1
                                ! fac = 0: qq1 returns the center position
      real, intent(out):: a1(2), a2(2), a3(2), a4(2)
! Local
      real qq1(3), qq2(3), qq3(3), qq4(3)
      real p1(3), p2(3), p3(3), p4(3)
      real ec(3)
      real(f_p):: dd, d1, d2, d3, d4
      integer k

! Transform to (x,y,z)
      call latlon2xyz(q1, p1)
      call latlon2xyz(q2, p2)
      call latlon2xyz(q3, p3)
      call latlon2xyz(q4, p4)

! Get center position:
      do k=1,3
         ec(k) = p1(k) + p2(k) + p3(k) + p4(k)
      enddo
      dd = sqrt( ec(1)**2 + ec(2)**2 + ec(3)**2 )

      do k=1,3
         ec(k) = ec(k) / dd   ! cell center position
      enddo

! Perform the "extrapolation" in 3D (x-y-z) 
      do k=1,3
         qq1(k) = ec(k) + fac*(p1(k)-ec(k)) 
         qq2(k) = ec(k) + fac*(p2(k)-ec(k)) 
         qq3(k) = ec(k) + fac*(p3(k)-ec(k)) 
         qq4(k) = ec(k) + fac*(p4(k)-ec(k)) 
      enddo

!--------------------------------------------------------
! Force the points to be on the sphere via normalization
!--------------------------------------------------------
      d1 = sqrt( qq1(1)**2 + qq1(2)**2 + qq1(3)**2 )
      d2 = sqrt( qq2(1)**2 + qq2(2)**2 + qq2(3)**2 )
      d3 = sqrt( qq3(1)**2 + qq3(2)**2 + qq3(3)**2 )
      d4 = sqrt( qq4(1)**2 + qq4(2)**2 + qq4(3)**2 )
      do k=1,3
         qq1(k) = qq1(k) / d1
         qq2(k) = qq2(k) / d2
         qq3(k) = qq3(k) / d3
         qq4(k) = qq4(k) / d4
      enddo

!----------------------------------------
! Transform back to lat-lon coordinates:
!----------------------------------------

      call cart_to_latlon(1, qq1, a1(1), a1(2))
      call cart_to_latlon(1, qq2, a2(1), a2(2))
      call cart_to_latlon(1, qq3, a3(1), a3(2))
      call cart_to_latlon(1, qq4, a4(1), a4(2))

 end subroutine expand_cell


 subroutine cell_center2(q1, q2, q3, q4, e2)
      real , intent(in ) :: q1(2), q2(2), q3(2), q4(2)
      real , intent(out) :: e2(2)
! Local
      real p1(3), p2(3), p3(3), p4(3)
      real ec(3)
      real dd
      integer k

      call latlon2xyz(q1, p1)
      call latlon2xyz(q2, p2)
      call latlon2xyz(q3, p3)
      call latlon2xyz(q4, p4)

      do k=1,3
         ec(k) = p1(k) + p2(k) + p3(k) + p4(k)
      enddo
      dd = sqrt( ec(1)**2 + ec(2)**2 + ec(3)**2 )

      do k=1,3
         ec(k) = ec(k) / dd
      enddo

      call cart_to_latlon(1, ec, e2(1), e2(2))

 end subroutine cell_center2


 subroutine cell_center3(p1, p2, p3, p4, ec)
! Get center position of a cell
         real , intent(IN)  :: p1(3), p2(3), p3(3), p4(3)
         real , intent(OUT) :: ec(3)
! Local
         real dd
         integer k

         do k=1,3
            ec(k) = p1(k) + p2(k) + p3(k) + p4(k)
         enddo
         dd = sqrt( ec(1)**2 + ec(2)**2 + ec(3)**2 )

         do k=1,3
            ec(k) = ec(k) / dd
         enddo

 end subroutine cell_center3



 real function get_area(p1, p4, p2, p3, radius)
!-----------------------------------------------
 real, intent(in), dimension(2):: p1, p2, p3, p4
 real, intent(in), optional:: radius
!-----------------------------------------------
 real e1(3), e2(3), e3(3)
 real ang1, ang2, ang3, ang4
 real pi

 pi = 4.*atan(1.)

! S-W: 1
       call latlon2xyz(p1, e1)   ! p1
       call latlon2xyz(p2, e2)   ! p2
       call latlon2xyz(p4, e3)   ! p4
       ang1 = spherical_angle(e1, e2, e3)
!----
! S-E: 2
!----
       call latlon2xyz(p2, e1)
       call latlon2xyz(p3, e2)
       call latlon2xyz(p1, e3)
       ang2 = spherical_angle(e1, e2, e3)
!----
! N-E: 3
!----
       call latlon2xyz(p3, e1)
       call latlon2xyz(p4, e2)
       call latlon2xyz(p2, e3)
       ang3 = spherical_angle(e1, e2, e3)
!----
! N-W: 4
!----
       call latlon2xyz(p4, e1)
       call latlon2xyz(p3, e2)
       call latlon2xyz(p1, e3)
       ang4 = spherical_angle(e1, e2, e3)

       if ( present(radius) ) then
            get_area = (ang1 + ang2 + ang3 + ang4 - 2.*pi) * radius**2
       else
            get_area = ang1 + ang2 + ang3 + ang4 - 2.*pi
       endif

 end function get_area



 real function spherical_angle(p1, p2, p3)
 
!           p3
!         /
!        /
!       p1 ---> angle
!         \
!          \
!           p2

 real p1(3), p2(3), p3(3)

 real (f_p):: e1(3), e2(3), e3(3)
 real (f_p):: px, py, pz
 real (f_p):: qx, qy, qz
 real (f_p):: angle, ddd
 integer n

  do n=1,3
     e1(n) = p1(n)
     e2(n) = p2(n)
     e3(n) = p3(n)
  enddo

!-------------------------------------------------------------------
! Page 41, Silverman's book on Vector Algebra; spherical trigonmetry
!-------------------------------------------------------------------
! Vector P:
   px = e1(2)*e2(3) - e1(3)*e2(2) 
   py = e1(3)*e2(1) - e1(1)*e2(3) 
   pz = e1(1)*e2(2) - e1(2)*e2(1) 
! Vector Q:
   qx = e1(2)*e3(3) - e1(3)*e3(2) 
   qy = e1(3)*e3(1) - e1(1)*e3(3) 
   qz = e1(1)*e3(2) - e1(2)*e3(1) 

   ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz)

   if ( ddd <= 0.0 ) then
        angle = 0.
   else
        ddd = (px*qx+py*qy+pz*qz) / sqrt(ddd)
        if ( abs(ddd)>1.) then
             angle = 2.*atan(1.0)    ! 0.5*pi
        else
             angle = acos( ddd )
        endif
   endif

   spherical_angle = angle

 end function spherical_angle


 real function cos_angle(p1, p2, p3)
! As spherical_angle, but returns the cos(angle)
 real p1(3), p2(3), p3(3)

 real (f_p):: e1(3), e2(3), e3(3)
 real (f_p):: px, py, pz
 real (f_p):: qx, qy, qz
 real (f_p):: angle, ddd
 integer n

  do n=1,3
     e1(n) = p1(n)
     e2(n) = p2(n)
     e3(n) = p3(n)
  enddo

!-------------------------------------------------------------------
! Page 41, Silverman's book on Vector Algebra; spherical trigonmetry
!-------------------------------------------------------------------
! Vector P:
   px = e1(2)*e2(3) - e1(3)*e2(2) 
   py = e1(3)*e2(1) - e1(1)*e2(3) 
   pz = e1(1)*e2(2) - e1(2)*e2(1) 
! Vector Q:
   qx = e1(2)*e3(3) - e1(3)*e3(2) 
   qy = e1(3)*e3(1) - e1(1)*e3(3) 
   qz = e1(1)*e3(2) - e1(2)*e3(1) 

   ddd = sqrt( (px**2+py**2+pz**2)*(qx**2+qy**2+qz**2) )

   if ( ddd > 0. ) then
        angle = (px*qx+py*qy+pz*qz) / ddd 
   else
        angle = 1.
   endif

   cos_angle = angle

 end function cos_angle



 real function g_sum(p, ifirst, ilast, jfirst, jlast, ngc, area, mode, reproduce)
! Fast version of globalsum 
      integer, intent(IN) :: ifirst, ilast
      integer, intent(IN) :: jfirst, jlast, ngc
      integer, intent(IN) :: mode  ! if ==1 divided by area
      logical, intent(in), optional :: reproduce
      real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
      real, intent(IN) :: area(ifirst-ngc:ilast+ngc,jfirst-ngc:jlast+ngc)
      integer :: i,j
      real gsum
         
      if ( .not. g_sum_initialized ) then
         global_area = mpp_global_sum(domain, area, flags=BITWISE_EXACT_SUM)
         if ( gid==0 ) write(*,*) 'Global Area=',global_area
         g_sum_initialized = .true.
      end if
 
!-------------------------
! FMS global sum algorithm:
!-------------------------
      if ( present(reproduce) ) then
         if (reproduce) then
            gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast), &
                                  flags=BITWISE_EXACT_SUM)
         else
            gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast))
         endif
      else
!-------------------------
! Quick local sum algorithm
!-------------------------
         gsum = 0.
         do j=jfirst,jlast
            do i=ifirst,ilast
               gsum = gsum + p(i,j)*area(i,j)
            enddo
         enddo
         call mp_reduce_sum(gsum)
      endif

      if ( mode==1 ) then
           g_sum = gsum / global_area
      else
           g_sum = gsum
      endif

 end function g_sum


 real function global_qsum(p, ifirst, ilast, jfirst, jlast)
! quick global sum without area weighting
      integer, intent(IN) :: ifirst, ilast
      integer, intent(IN) :: jfirst, jlast
      real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
      integer :: i,j
      real gsum
         
      gsum = 0.
      do j=jfirst,jlast
         do i=ifirst,ilast
            gsum = gsum + p(i,j)
         enddo
      enddo
      call mp_reduce_sum(gsum)

      global_qsum  = gsum

 end function global_qsum


 subroutine global_mx(q, n_g, qmin, qmax)
     integer, intent(in):: n_g
     real, intent(in)::q(is-n_g:ie+n_g, js-n_g:je+n_g)
     real, intent(out):: qmin, qmax
     integer i,j

      qmin = q(is,js)
      qmax = qmin
      do j=js,je
         do i=is,ie
            qmin = min(qmin, q(i,j))
            qmax = max(qmax, q(i,j))
         enddo
      enddo
      call mp_reduce_min(qmin)
      call mp_reduce_max(qmax)

 end subroutine global_mx

 subroutine global_mx_c(q, i1, i2, j1, j2, qmin, qmax)
! For computing global max/min at cell Corners
     integer, intent(in):: i1, i2, j1, j2
     real, intent(in)   :: q(i1:i2,j1:j2)
     real, intent(out)  :: qmin, qmax
     integer i,j

      qmin = q(i1,j1)
      qmax = qmin
      do j=j1,j2
         do i=i1,i2
            qmin = min(qmin, q(i,j))
            qmax = max(qmax, q(i,j))
         enddo
      enddo
      call mp_reduce_min(qmin)
      call mp_reduce_max(qmax)

 end subroutine global_mx_c



  subroutine fill_ghost(q, npx, npy, value)
  real, intent(inout):: q(isd:ied,jsd:jed)
  integer, intent(in):: npx, npy
  real, intent(in):: value
  integer i,j

     do j=jsd,jed
        do i=isd,ied
           if ( (i<1 .and. j<1) ) then
                q(i,j) = value
           endif
           if ( i>(npx-1) .and. j<1 ) then
                q(i,j) = value
           endif
           if ( i>(npx-1) .and. j>(npy-1) ) then
                q(i,j) = value
           endif
           if ( i<1 .and. j>(npy-1) ) then
                q(i,j) = value
           endif
        enddo
     enddo

  end subroutine fill_ghost



 subroutine make_eta_level(km, pe, area, kks, ak, bk)
  integer, intent(in ):: km
  integer, intent(out):: kks
  real, intent(in):: area(isd:ied,jsd:jed)
  real, intent(in):: pe(is-1:ie+1,km+1,js-1:je+1)
  real, intent(out):: ak(km+1), bk(km+1)
! local:
  real ph(km+1)
  real, allocatable:: pem(:,:)
  real(kind=4) :: p4
  integer k, i, j

     ph(1) = pe(is,1,js)

     allocate ( pem(is:ie,js:je) )

! Compute global mean values:
     do k=2,km+1
        do j=js,je
           do i=is,ie
               pem(i,j) = pe(i,k,j)
           enddo
        enddo
! Make it the same across all PEs
!       p4 = g_sum(pem, is, ie, js, je, ng, area, mode=1)
        p4 = g_sum(pem, is, ie, js, je, ng, area, 1)
        ph(k) = p4
     enddo

! Faking a's and b's for code compatibility with hybrid sigma-p
     kks = 0
     ak(1) = ph(1)
     bk(1) = 0.
     ak(km+1) = 0.
     bk(km+1) = 1.

     do k=2,km
        bk(k) = (ph(k) - ph(1)) / (ph(km+1)-ph(1))
        ak(k) = ph(1)*(1.-bk(k))
     enddo

    if ( gid==0 ) then
         write(*,*) 'Make_eta_level ...., ptop=', ph(1)
         ptop = ph(1)
#ifdef PRINT_GRID
         do k=1,km+1
            write(*,*) ph(k), ak(k), bk(k)
         enddo
#endif
    endif

    deallocate ( pem )

 end subroutine make_eta_level

 subroutine invert_matrix(n, a, x)
  integer, intent (in) :: n
  integer :: i,j,k
  real, intent (inout), dimension (n,n):: a
  real, intent (out), dimension (n,n):: x   ! inverted maxtrix
  real, dimension (n,n) :: b
  integer indx(n)
 
  do i = 1, n
     do j = 1, n
        b(i,j) = 0.0
     end do
  end do

  do i = 1, n
     b(i,i) = 1.0
  end do
 
  call elgs (a,n,indx)
 
  do i = 1, n-1
     do j = i+1, n
        do k = 1, n
           b(indx(j),k) = b(indx(j),k) - a(indx(j),i)*b(indx(i),k)
        end do
     end do
  end do
 
  do i = 1, n
     x(n,i) = b(indx(n),i)/a(indx(n),n)
     do j = n-1, 1, -1
        x(j,i) = b(indx(j),i)
        do k = j+1, n
           x(j,i) = x(j,i)-a(indx(j),k)*x(k,i)
        end do
        x(j,i) =  x(j,i)/a(indx(j),j)
     end do
  end do

 end subroutine invert_matrix
 

 subroutine elgs (a,n,indx)

!------------------------------------------------------------------
! subroutine to perform the partial-pivoting gaussian elimination.
! a(n,n) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
!------------------------------------------------------------------
 
  integer, intent (in) :: n
  integer :: i,j,k,itmp
  integer, intent (out), dimension (n) :: indx
  real, intent (inout), dimension (n,n) :: a
!
  real :: c1,pi,pi1,pj
  real, dimension (n) :: c
 
  do i = 1, n
     indx(i) = i
  end do
!
! find the rescaling factors, one from each row
!
  do i = 1, n
     c1= 0.0
     do j = 1, n
        c1 = max(c1,abs(a(i,j)))
     end do
     c(i) = c1
  end do
!
! search the pivoting (largest) element from each column
!
  do j = 1, n-1
     pi1 = 0.0
     do i = j, n
        pi = abs(a(indx(i),j))/c(indx(i))
        if (pi > pi1) then
            pi1 = pi
            k   = i
        endif
     end do
!
! interchange the rows via indx(n) to record pivoting order
!
    itmp    = indx(j)
    indx(j) = indx(k)
    indx(k) = itmp
    do i = j+1, n
       pj  = a(indx(i),j)/a(indx(j),j)
!
! record pivoting ratios below the diagonal
!
       a(indx(i),j) = pj
!
! modify other elements accordingly
!
       do k = j+1, n
          a(indx(i),k) = a(indx(i),k)-pj*a(indx(j),k)
       end do
     end do
  end do
 
 end subroutine elgs

 end module fv_grid_utils_mod

