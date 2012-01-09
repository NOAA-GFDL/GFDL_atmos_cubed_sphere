 module test_cases_mod

      use constants_mod,     only: radius, pi, omega, grav, kappa, rdgas, cp_air
      use init_hydro_mod,    only: p_var, hydro_eq
      use fv_mp_mod,         only: gid, masterproc, domain, tile, ng,         &
                                   is,js,ie,je, isd,jsd,ied,jed, &
                                   domain_decomp, fill_corners, XDir, YDir, &
                                   mp_stop, mp_reduce_sum, mp_reduce_max, mp_gather, mp_bcst
      use fv_grid_utils_mod, only: cubed_to_latlon, great_circle_dist, mid_pt_sphere,   &
                                   ks, ptop, ptop_min, fC, f0, deglat, inner_prod, normalize_vect, &
                                   ee1, ee2, ew, es, g_sum, latlon2xyz, cart_to_latlon, make_eta_level
      use fv_surf_map_mod,   only: surfdrv

      use fv_grid_tools_mod, only: grid, agrid, cubed_sphere, latlon,  todeg, missing,  &
                                   dx,dy, dxa,dya, rdxa, rdya, dxc,dyc, area, rarea,rarea_c, &
                                   ctoa, atod, dtoa, atoc, atob_s, mp_update_dwinds, rotate_winds, &
                                   globalsum, get_unit_vector, unit_vect2,                         &
                                   dx_const, dy_const
      use fv_eta_mod,        only: compute_dz_L32, compute_dz_L101, set_hybrid_z, gw_1d

      use mpp_mod,           only: mpp_error, FATAL
      use mpp_domains_mod,   only: mpp_update_domains
      use mpp_parameter_mod, only: AGRID_PARAM=>AGRID,CGRID_NE_PARAM=>CGRID_NE, &
                                   SCALAR_PAIR
      use fv_sg_mod,         only: qsmith
!     use fv_diagnostics_mod, only: prt_maxmin


      implicit none
      private

! Test Case Number  
!                   -1 = Divergence conservation test
!                    0 = Idealized non-linear deformational flow
!                    1 = Cosine Bell advection
!                    2 = Zonal geostrophically balanced flow
!                    3 = non-rotating potential flow 
!                    4 = Tropical cyclones (merger of Rankine vortices)
!                    5 = Zonal geostrophically balanced flow over an isolated mountain
!                    6 = Rossby Wave number 4 
!                    7 = Barotropic instability
!                    8 = Potential flow (as in 5 but no rotation and initially at rest)
!                    9 = Polar vortex
!                   10 = hydrostatically balanced 3D test with idealized mountain
!                   11 = Use this for cold starting the climate model with USGS terrain
!                   12 = Jablonowski & Williamson Baroclinic test case (Steady State)
!                   13 = Jablonowski & Williamson Baroclinic test case Perturbation
!                   14 = Use this for cold starting the Aqua-planet model
!                   15 = Small Earth density current
!                   16 = 3D hydrostatic non-rotating Gravity waves
!                   17 = 3D hydrostatic rotating Inertial Gravity waves (case 6-3-0)
!                   18 = 3D mountain-induced Rossby wave
!                   19 = As in 15 but without rotation
!                   20 = 3D non-hydrostatic lee vortices; non-rotating (small planet)
!                   21 = 3D non-hydrostatic lee vortices; rotating     (small planet)
!                  101 = 3D non-hydrostatic Large-Eddy-Simulation (LES) with hybrid_z IC
      integer :: test_case
! alpha = angle of axis rotation about the poles
      real   :: alpha = 0.
! Ubar = initial wind speed parameter
      real   :: Ubar
! gh0 = initial surface height parameter
      real   :: gh0

! Case 0 parameters
      real :: p0 = 3.0
      real :: rgamma = 5.0
      real :: lat0 = pi/2.0 !pi/4.8
      real :: lon0 = 0.0 !pi-0.8

!  pi_shift moves the initial location of the cosine bell for Case 1
      real, parameter :: pi_shift = 0.0 !3.0*pi/4.

!  case 9 parameters 
      real  , allocatable :: case9_B(:,:)
      real   :: AofT(2)

!  Validating fields used in statistics
      real  , allocatable :: phi0(:,:,:) ! Validating Field
      real  , allocatable :: ua0(:,:,:)  ! Validating U-Wind
      real  , allocatable :: va0(:,:,:)  ! Validating V-Wind
      real  , allocatable :: gh_table(:), lats_table(:)

!  Initial Conservation statistics ; total mass ; enstrophy ; energy
      real   :: tmass_orig
      real   :: tvort_orig
      real   :: tener_orig

 ! -1:null_op, 0:All-Grids, 1:C-Grid, 2:D-Grid, 3:A-Grid, 4:A-Grid then Rotate, 5:D-Grid with unit vectors then Rotate
      integer, parameter :: initWindsCase0 =-1 
      integer, parameter :: initWindsCase1 = 1
      integer, parameter :: initWindsCase2 = 5 
      integer, parameter :: initWindsCase5 = 5
      integer, parameter :: initWindsCase6 =-1 
      integer, parameter :: initWindsCase9 =-1

      public :: test_case, alpha
      public :: init_case, get_stats, check_courant_numbers, output, output_ncdf
      public :: case9_forcing1, case9_forcing2
      public :: init_double_periodic, init_latlon

 !---- version number -----
      character(len=128) :: version = '$Id: test_cases.F90,v 19.0 2012/01/06 19:59:28 fms Exp $'
      character(len=128) :: tagname = '$Name: siena $'

      contains

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     init_winds :: initialize the winds 
!
      subroutine init_winds(UBar, u,v,ua,va,uc,vc, defOnGrid, npx, npy, ng, ndims, nregions)
 ! defOnGrid = -1:null_op, 0:All-Grids, 1:C-Grid, 2:D-Grid, 3:A-Grid, 4:A-Grid then Rotate, 5:D-Grid with unit vectors then Rotate

      real  ,    intent(INOUT) :: UBar
      real ,      intent(INOUT) ::    u(isd:ied  ,jsd:jed+1)
      real ,      intent(INOUT) ::    v(isd:ied+1,jsd:jed  )
      real ,      intent(INOUT) ::   uc(isd:ied+1,jsd:jed  )
      real ,      intent(INOUT) ::   vc(isd:ied  ,jsd:jed+1)
      real ,      intent(INOUT) ::   ua(isd:ied  ,jsd:jed  )
      real ,      intent(INOUT) ::   va(isd:ied  ,jsd:jed  )
      integer,      intent(IN) :: defOnGrid
      integer,      intent(IN) :: npx, npy
      integer,      intent(IN) :: ng
      integer,      intent(IN) :: ndims
      integer,      intent(IN) :: nregions

      real   :: p1(2),p2(2),p3(2),p4(2), pt(2)
      real :: e1(3), e2(3), ex(3), ey(3)

      real   :: dist, r, r0 
      integer :: i,j,k,n
      real :: utmp, vtmp

      real :: psi_b(isd:ied+1,jsd:jed+1), psi(isd:ied,jsd:jed), psi1, psi2 

 200  format(i4.4,'x',i4.4,'x',i4.4,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14)

      psi(:,:) = 1.e25
      psi_b(:,:) = 1.e25
      do j=jsd,jed
         do i=isd,ied
            psi(i,j) = (-1.0 * Ubar * radius *( sin(agrid(i,j,2))                  *cos(alpha) - &
                                            cos(agrid(i,j,1))*cos(agrid(i,j,2))*sin(alpha) ) )
         enddo
      enddo
      call mpp_update_domains( psi, domain )
      do j=jsd,jed+1
         do i=isd,ied+1
            psi_b(i,j) = (-1.0 * Ubar * radius *( sin(grid(i,j,2))                 *cos(alpha) - &
                                              cos(grid(i,j,1))*cos(grid(i,j,2))*sin(alpha) ) )
         enddo
      enddo

      if ( (cubed_sphere) .and. (defOnGrid==0) ) then
         do j=js,je+1
            do i=is,ie
               dist = dx(i,j)
               vc(i,j) = (psi_b(i+1,j)-psi_b(i,j))/dist
               if (dist==0) vc(i,j) = 0.
            enddo
         enddo
         do j=js,je
            do i=is,ie+1
               dist = dy(i,j)
               uc(i,j) = -1.0*(psi_b(i,j+1)-psi_b(i,j))/dist
               if (dist==0) uc(i,j) = 0.
            enddo
         enddo
         call mpp_update_domains( uc, vc, domain, gridtype=CGRID_NE_PARAM)
         call fill_corners(uc, vc, npx, npy, VECTOR=.true., CGRID=.true.)
         do j=js,je
            do i=is,ie+1
               dist = dxc(i,j)
               v(i,j) = (psi(i,j)-psi(i-1,j))/dist
               if (dist==0) v(i,j) = 0.
            enddo
         enddo
         do j=js,je+1
            do i=is,ie
               dist = dyc(i,j)
               u(i,j) = -1.0*(psi(i,j)-psi(i,j-1))/dist
               if (dist==0) u(i,j) = 0.
            enddo
         enddo
         call mp_update_dwinds(u, v, npx, npy)
         do j=js,je
            do i=is,ie
               psi1 = 0.5*(psi(i,j)+psi(i,j-1))
               psi2 = 0.5*(psi(i,j)+psi(i,j+1))
               dist = dya(i,j)
               ua(i,j) = -1.0 * (psi2 - psi1) / (dist)
               if (dist==0) ua(i,j) = 0.
               psi1 = 0.5*(psi(i,j)+psi(i-1,j))
               psi2 = 0.5*(psi(i,j)+psi(i+1,j))
               dist = dxa(i,j)
               va(i,j) = (psi2 - psi1) / (dist)
               if (dist==0) va(i,j) = 0.
            enddo
         enddo

      elseif ( (cubed_sphere) .and. (defOnGrid==1) ) then
         do j=js,je+1
            do i=is,ie
               dist = dx(i,j)
               vc(i,j) = (psi_b(i+1,j)-psi_b(i,j))/dist
               if (dist==0) vc(i,j) = 0.
            enddo
         enddo
         do j=js,je
            do i=is,ie+1
               dist = dy(i,j)
               uc(i,j) = -1.0*(psi_b(i,j+1)-psi_b(i,j))/dist
               if (dist==0) uc(i,j) = 0.
            enddo
         enddo
         call mpp_update_domains( uc, vc, domain, gridtype=CGRID_NE_PARAM)
         call fill_corners(uc, vc, npx, npy, VECTOR=.true., CGRID=.true.)
         call ctoa(uc,vc,ua,va,npx,npy,ng)
         call atod(ua,va,u ,v ,npx,npy,ng)
        ! call d2a2c(npx,npy,1, is,ie, js,je, ng, u(isd,jsd),v(isd,jsd), &
        !            ua(isd,jsd),va(isd,jsd), uc(isd,jsd),vc(isd,jsd))
      elseif ( (cubed_sphere) .and. (defOnGrid==2) ) then
         do j=js,je
            do i=is,ie+1
               dist = dxc(i,j)
               v(i,j) = (psi(i,j)-psi(i-1,j))/dist
               if (dist==0) v(i,j) = 0.            
            enddo
         enddo
         do j=js,je+1
            do i=is,ie
               dist = dyc(i,j)
               u(i,j) = -1.0*(psi(i,j)-psi(i,j-1))/dist
               if (dist==0) u(i,j) = 0. 
            enddo
         enddo
         call mp_update_dwinds(u, v, npx, npy)
         call dtoa( u, v,ua,va,npx,npy,ng)
         call atoc(ua,va,uc,vc,npx,npy,ng) 
      elseif ( (cubed_sphere) .and. (defOnGrid==3) ) then
         do j=js,je
            do i=is,ie
               psi1 = 0.5*(psi(i,j)+psi(i,j-1))
               psi2 = 0.5*(psi(i,j)+psi(i,j+1))
               dist = dya(i,j)
               ua(i,j) = -1.0 * (psi2 - psi1) / (dist)
               if (dist==0) ua(i,j) = 0.
               psi1 = 0.5*(psi(i,j)+psi(i-1,j))
               psi2 = 0.5*(psi(i,j)+psi(i+1,j))
               dist = dxa(i,j)
               va(i,j) = (psi2 - psi1) / (dist)
               if (dist==0) va(i,j) = 0.
            enddo
         enddo
         call mpp_update_domains( ua, va, domain, gridtype=AGRID_PARAM)
         call atod(ua,va, u, v,npx,npy,ng)
         call atoc(ua,va,uc,vc,npx,npy,ng)
      elseif ( (latlon) .or. (defOnGrid==4) ) then

         do j=js,je
            do i=is,ie
               ua(i,j) =  Ubar * ( COS(agrid(i,j,2))*COS(alpha) + &
                                     SIN(agrid(i,j,2))*COS(agrid(i,j,1))*SIN(alpha) )
               va(i,j) = -Ubar *   SIN(agrid(i,j,1))*SIN(alpha)  
               call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p1)
               call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p2)
               call mid_pt_sphere(grid(i+1,j,1:2), grid(i+1,j+1,1:2), p3)
               call mid_pt_sphere(grid(i,j+1,1:2), grid(i+1,j+1,1:2), p4)
               if (cubed_sphere) call rotate_winds(ua(i,j), va(i,j), p1,p2,p3,p4, agrid(i,j,1:2), 2, 1)

               psi1 = 0.5*(psi(i,j)+psi(i,j-1))
               psi2 = 0.5*(psi(i,j)+psi(i,j+1))
               dist = dya(i,j)
    if ( (tile==1) .and.(i==1) ) print*, ua(i,j), -1.0 * (psi2 - psi1) / (dist)

            enddo
         enddo
         call mpp_update_domains( ua, va, domain, gridtype=AGRID_PARAM)
         call atod(ua,va, u, v,npx,npy,ng)
         call atoc(ua,va,uc,vc,npx,npy,ng)

     elseif ( (latlon) .or. (defOnGrid==5) ) then
! SJL mods:
! v-wind:
         do j=js,je
            do i=is,ie+1
               p1(:) = grid(i  ,j ,1:2)
               p2(:) = grid(i,j+1 ,1:2)
               call mid_pt_sphere(p1, p2, pt)
               call get_unit_vector(p1, p2, e2)
!              call unit_vect2(p1, p2, e2)
               call get_latlon_vector(pt, ex, ey)
               utmp =  Ubar * ( COS(pt(2))*COS(alpha) + &
                                SIN(pt(2))*COS(pt(1))*SIN(alpha) )
               vtmp = -Ubar *   SIN(pt(1))*SIN(alpha)
               v(i,j) = utmp*inner_prod(e2,ex) + vtmp*inner_prod(e2,ey)
            enddo
         enddo
! D grid u-wind:
         do j=js,je+1
            do i=is,ie
               p1(:) = grid(i  ,j  ,1:2)
               p2(:) = grid(i+1,j  ,1:2)
               call mid_pt_sphere(p1, p2, pt)
               call get_unit_vector(p1, p2, e1)
!              call unit_vect2(p1, p2, e1)
               call get_latlon_vector(pt, ex, ey)
               utmp =  Ubar * ( COS(pt(2))*COS(alpha) + &
                                SIN(pt(2))*COS(pt(1))*SIN(alpha) )
               vtmp = -Ubar *   SIN(pt(1))*SIN(alpha)
               u(i,j) = utmp*inner_prod(e1,ex) + vtmp*inner_prod(e1,ey)
            enddo
         enddo

         call mp_update_dwinds(u, v, npx, npy)
         call dtoa( u, v,ua,va,npx,npy,ng)
         call atoc(ua,va,uc,vc,npx,npy,ng)
     else
         !print*, 'Choose an appropriate grid to define the winds on'
         !stop
     endif

      end subroutine init_winds
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     init_case :: initialize the Williamson test cases:
!                  case 1 (2-D advection of a cosine bell)
!                  case 2 (Steady State Zonal Geostrophic Flow)
!                  case 5 (Steady State Zonal Geostrophic Flow over Mountain)
!                  case 6 (Rossby Wave-4 Case)
!                  case 9 (Stratospheric Vortex Breaking Case)
!
      subroutine init_case(u,v,w,pt,delp,q,phis, ps,pe,peln,pk,pkz,  uc,vc, ua,va, ak, bk,  &
                           npx, npy, npz, ng, ncnst, nwat, k_top, ndims, nregions,        &
                           dry_mass, mountain, moist_phys, hydrostatic, hybrid_z, delz, ze0)

      real ,      intent(INOUT) ::    u(isd:ied  ,jsd:jed+1,npz)
      real ,      intent(INOUT) ::    v(isd:ied+1,jsd:jed  ,npz)
      real ,      intent(INOUT) ::    w(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(INOUT) ::   pt(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(INOUT) :: delp(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(INOUT) ::    q(isd:ied  ,jsd:jed  ,npz, ncnst)

      real ,      intent(INOUT) :: phis(isd:ied  ,jsd:jed  )

      real ,      intent(INOUT) ::   ps(isd:ied  ,jsd:jed  )
      real ,      intent(INOUT) ::   pe(is-1:ie+1,npz+1,js-1:je+1)
      real ,      intent(INOUT) ::   pk(is:ie    ,js:je    ,npz+1)
      real ,      intent(INOUT) :: peln(is :ie   ,npz+1    ,js:je)
      real ,      intent(INOUT) ::  pkz(is:ie    ,js:je    ,npz  )

      real ,      intent(INOUT) ::   uc(isd:ied+1,jsd:jed  ,npz)
      real ,      intent(INOUT) ::   vc(isd:ied  ,jsd:jed+1,npz)
      real ,      intent(INOUT) ::   ua(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(INOUT) ::   va(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(inout) :: delz(is:ie,js:je,npz)
      real ,      intent(inout)   ::  ze0(is:ie,js:je,npz+1)

      real ,      intent(inout) ::   ak(npz+1)
      real ,      intent(inout) ::   bk(npz+1)

      integer,      intent(IN) :: npx, npy, npz
      integer,      intent(IN) :: ng, ncnst, nwat
      integer,      intent(IN) :: k_top
      integer,      intent(IN) :: ndims
      integer,      intent(IN) :: nregions

      real,         intent(IN) :: dry_mass
      logical,      intent(IN) :: mountain
      logical,      intent(IN) :: moist_phys
      logical,      intent(IN) :: hydrostatic
      logical,      intent(IN) :: hybrid_z

      real   ::  tmp(1-ng:npx  +ng,1-ng:npy  +ng,1:nregions)
      real   :: tmp1(1   :npx     ,1   :npy     ,1:nregions)

      real   :: p1(2)      ! Temporary Point
      real   :: p2(2)      ! Temporary Point
      real   :: p3(2)      ! Temporary Point
      real   :: p4(2)      ! Temporary Point
      real   :: pa(2)      ! Temporary Point
      real   :: pb(2)      ! Temporary Point
      real   :: pcen(2)    ! Temporary Point
      real   :: e1(3), e2(3), e3(3), ex(3), ey(3)
      real   :: dist, r, r0, omg, A, B, C
      integer :: i,j,k,nreg,z,zz
      integer :: i0,j0,n0
      real   :: utmp,vtmp,ftmp
      real   :: rk

      integer, parameter :: jm = 5761
      real   :: ll_phi(jm)
      real   ::   ll_u(jm)
      real   ::   ll_j(jm)
      real   ::   cose(jm)
      real   ::   sine(jm)
      real   ::   cosp(jm)
      real   :: ddeg, deg, DDP, DP, ph5
      real   :: myB, myC, yy
      integer   :: jj,jm1

      real :: Vtx, p, w_p
      real :: x1,y1,z1,x2,y2,z2,ang

      integer :: initWindsCase

      real :: dummy
      real :: ftop
      real :: v1,v2
      real :: m=1
      real :: n=1
      real :: L1_norm
      real :: L2_norm
      real :: Linf_norm
      real :: pmin, pmin1
      real :: pmax, pmax1
      real :: grad(isd:ied  ,jsd:jed,2)
      real :: div0(isd:ied  ,jsd:jed  ) 
      real :: vor0(isd:ied  ,jsd:jed  )
      real :: divg(isd:ied  ,jsd:jed  )
      real :: vort(isd:ied  ,jsd:jed  )
      real :: ztop, rgrav, p00, pturb, zmid, pk0, t00
      real :: dz1(npz), ppt(npz)
      real :: ze1(npz+1), pe1(npz+1)

      integer :: nlon,nlat
      character(len=80) :: oflnm, hgtflnm

! Baroclinic Test Case 12
      real :: eta, eta_0, eta_s, eta_t
      real :: eta_v(npz), press, anti_rot
      real :: T_0, T_mean, delta_T, lapse_rate, n2, zeta, s0
      real :: pt1,pt2,pt3,pt4,pt5,pt6, pt7, pt8, pt9, u1, pt0
      real :: uu1, uu2, uu3, vv1, vv2, vv3
!     real wbuffer(npx+1,npz)
!     real sbuffer(npy+1,npz)
      real wbuffer(npy+2,npz)
      real sbuffer(npx+2,npz)

      allocate ( phi0(isd:ied  ,jsd:jed  ,npz) )
      allocate (  ua0(isd:ied  ,jsd:jed  ,npz) )
      allocate (  va0(isd:ied  ,jsd:jed  ,npz) )

      pe(:,:,:) = 0.0
      pt(:,:,:) = 1.0
      f0(:,:) = huge(dummy)
      fC(:,:) = huge(dummy)
      do j=jsd,jed+1
         do i=isd,ied+1
            fC(i,j) = 2.*omega*( -1.*cos(grid(i,j,1))*cos(grid(i,j,2))*sin(alpha) + &
                                     sin(grid(i,j,2))*cos(alpha) )
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            f0(i,j) = 2.*omega*( -1.*cos(agrid(i,j,1))*cos(agrid(i,j,2))*sin(alpha) + &
                                     sin(agrid(i,j,2))*cos(alpha) )
         enddo
      enddo
      call mpp_update_domains( f0, domain )
      if (cubed_sphere) call fill_corners(f0, npx, npy, YDir)

      delp(isd:is-1,jsd:js-1,1:npz)=0.
      delp(isd:is-1,je+1:jed,1:npz)=0.
      delp(ie+1:ied,jsd:js-1,1:npz)=0.
      delp(ie+1:ied,je+1:jed,1:npz)=0.

#if defined(SW_DYNAMICS)
      select case (test_case)
      case(-2)
      case(-1)
         Ubar = (2.0*pi*radius)/(12.0*86400.0)
         gh0  = 2.94e4
         phis = 0.0
         do j=js,je
            do i=is,ie
               delp(i,j,1) = gh0 - (radius*omega*Ubar + (Ubar*Ubar)/2.) * &
                             ( -1.*cos(agrid(i  ,j  ,1))*cos(agrid(i  ,j  ,2))*sin(alpha) + &
                                   sin(agrid(i  ,j  ,2))*cos(alpha) ) ** 2.0
            enddo
         enddo
         call init_winds(UBar, u,v,ua,va,uc,vc, 1, npx, npy, ng, ndims, nregions)

! Test Divergence operator at cell centers
         do j=js,je
            do i=is,ie
               divg(i,j) = (rarea(i,j)) * ( (uc(i+1,j,1)*dy(i+1,j) - uc(i,j,1)*dy(i,j)) + &
                                            (vc(i,j+1,1)*dx(i,j+1) - vc(i,j,1)*dx(i,j)) )
      if ( (tile==1) .and. (i==1) ) write(*,200) i,j,tile, divg(i,j), uc(i,j,1), uc(i+1,j,1), vc(i,j,1), vc(i,j+1,1)
            enddo
         enddo
! Test Vorticity operator at cell centers
         do j=js,je
            do i=is,ie
               vort(i,j) = (rarea(i,j)) * ( (v(i+1,j,1)*dy(i+1,j) - v(i,j,1)*dy(i,j)) - &
                                            (u(i,j+1,1)*dx(i,j+1) - u(i,j,1)*dx(i,j)) )
           enddo
        enddo
        div0(:,:) = 1.e-20
     ! call mpp_update_domains( div0, domain )
     ! call mpp_update_domains( vor0, domain )
     ! call mpp_update_domains( divg, domain )
     ! call mpp_update_domains( vort, domain )
      call get_scalar_stats( divg, div0, npx, npy, ndims, nregions, &
                             pmin, pmax, L1_norm, L2_norm, Linf_norm)
 200  format(i4.4,'x',i4.4,'x',i4.4,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14)
 201  format('          ',A,e21.14,' ',e21.14)
 202  format('          ',A,i4.4,'x',i4.4,'x',i4.4)
      if ( (gid == masterproc) ) then
          write(*,*) ' Error Norms of Analytical Divergence field C-Winds initialized'
          write(*,201) 'Divergence MAX error     : ', pmax
          write(*,201) 'Divergence MIN error     : ', pmin
          write(*,201) 'Divergence L1_norm       : ', L1_norm
          write(*,201) 'Divergence L2_norm       : ', L2_norm
          write(*,201) 'Divergence Linf_norm     : ', Linf_norm
      endif 

         call init_winds(UBar, u,v,ua,va,uc,vc, 3, npx, npy, ng, ndims, nregions)
! Test Divergence operator at cell centers
         do j=js,je
            do i=is,ie
               divg(i,j) = (rarea(i,j)) * ( (uc(i+1,j,1)*dy(i+1,j) - uc(i,j,1)*dy(i,j)) + &
                                            (vc(i,j+1,1)*dx(i,j+1) - vc(i,j,1)*dx(i,j)) )
      if ( (tile==1) .and. (i==1) ) write(*,200) i,j,tile, divg(i,j), uc(i,j,1), uc(i+1,j,1), vc(i,j,1), vc(i,j+1,1)
            enddo
         enddo
! Test Vorticity operator at cell centers
         do j=js,je
            do i=is,ie
               vort(i,j) = (rarea(i,j)) * ( (v(i+1,j,1)*dy(i+1,j) - v(i,j,1)*dy(i,j)) - &
                                            (u(i,j+1,1)*dx(i,j+1) - u(i,j,1)*dx(i,j)) )
           enddo
        enddo
        ua0 = ua
        va0 = va
        div0(:,:) = 1.e-20
      call get_scalar_stats( divg, div0, npx, npy, ndims, nregions, &
                             pmin, pmax, L1_norm, L2_norm, Linf_norm)
      if ( (gid == masterproc) ) then
          write(*,*) ' Error Norms of Analytical Divergence field A-Winds initialized'
          write(*,201) 'Divergence MAX error     : ', pmax
          write(*,201) 'Divergence MIN error     : ', pmin
          write(*,201) 'Divergence L1_norm       : ', L1_norm
          write(*,201) 'Divergence L2_norm       : ', L2_norm
          write(*,201) 'Divergence Linf_norm     : ', Linf_norm
      endif

         call init_winds(UBar, u,v,ua,va,uc,vc, 2, npx, npy, ng, ndims, nregions)
         !call d2a2c(npx,npy,1, is,ie, js,je, ng, u(isd,jsd,1),v(isd,jsd,1), &
         !           ua(isd,jsd,1),va(isd,jsd,1), uc(isd,jsd,1),vc(isd,jsd,1))
! Test Divergence operator at cell centers
         do j=js,je
            do i=is,ie
               divg(i,j) = (rarea(i,j)) * ( (uc(i+1,j,1)*dy(i+1,j) - uc(i,j,1)*dy(i,j)) + &
                                            (vc(i,j+1,1)*dx(i,j+1) - vc(i,j,1)*dx(i,j)) )
      if ( (tile==1) .and. ((i==1) .or.(i==npx-1)) ) write(*,200) i,j,tile, divg(i,j), uc(i,j,1), uc(i+1,j,1), vc(i,j,1), vc(i,j+1,1)
            enddo
         enddo
! Test Vorticity operator at cell centers
         do j=js,je
            do i=is,ie
               vort(i,j) = (rarea(i,j)) * ( (v(i+1,j,1)*dy(i+1,j) - v(i,j,1)*dy(i,j)) - &
                                            (u(i,j+1,1)*dx(i,j+1) - u(i,j,1)*dx(i,j)) )
           enddo
        enddo
        div0(:,:) = 1.e-20
      call get_scalar_stats( divg, div0, npx, npy, ndims, nregions, &
                             pmin, pmax, L1_norm, L2_norm, Linf_norm)
      if ( (gid == masterproc) ) then
          write(*,*) ' Error Norms of Analytical Divergence field D-Winds initialized'
          write(*,201) 'Divergence MAX error     : ', pmax
          write(*,201) 'Divergence MIN error     : ', pmin
          write(*,201) 'Divergence L1_norm       : ', L1_norm
          write(*,201) 'Divergence L2_norm       : ', L2_norm
          write(*,201) 'Divergence Linf_norm     : ', Linf_norm
      endif

      call mp_stop()
      stop
      case(0)
         do j=jsd,jed
            do i=isd,ied

               x1 = agrid(i,j,1) 
               y1 = agrid(i,j,2)
               z1 = radius

               p = p0 * cos(y1)
               Vtx = ((3.0*SQRT(2.0))/2.0) * (( 1.0/cosh(p) )**2.0) * tanh(p)
               w_p = 0.0
               if (p /= 0.0) w_p = Vtx/p 
               delp(i,j,1) = 1.0 - tanh( (p/rgamma) * sin(x1 - w_p*0.0) )
               ua(i,j,1) = w_p*(sin(lat0)*cos(agrid(i,j,2)) + cos(lat0)*cos(agrid(i,j,1) - lon0)*sin(agrid(i,j,2)))
               va(i,j,1) = w_p*cos(lat0)*sin(agrid(i,j,1) - lon0)
               ua(i,j,1) = ua(i,j,1)*radius/86400.0
               va(i,j,1) = va(i,j,1)*radius/86400.0

               call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p1)
               call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p2)
               call mid_pt_sphere(grid(i+1,j,1:2), grid(i+1,j+1,1:2), p3)
               call mid_pt_sphere(grid(i,j+1,1:2), grid(i+1,j+1,1:2), p4)      
               if (cubed_sphere) call rotate_winds(ua(i,j,1),va(i,j,1), p1,p2,p3,p4, agrid(i,j,1:2), 2, 1)

            enddo
         enddo
         call mpp_update_domains( ua, va, domain, gridtype=AGRID_PARAM)
         call atod(ua,va, u, v,npx,npy,ng)
         call mp_update_dwinds(u, v, npx, npy, npz)
         call atoc(ua,va,uc,vc,npx,npy,ng)
         call mpp_update_domains( uc, vc, domain, gridtype=CGRID_NE_PARAM)
         call fill_corners(uc, vc, npx, npy, npz, VECTOR=.true., CGRID=.true.)
         initWindsCase=initWindsCase0
      case(1)
         Ubar = (2.0*pi*radius)/(12.0*86400.0)
         gh0  = 1.0
         phis = 0.0
         r0 = radius/3. !RADIUS radius/3.
         p1(1) = pi/2. + pi_shift
         p1(2) = 0.
         do j=jsd,jed
            do i=isd,ied
               p2(1) = agrid(i,j,1)
               p2(2) = agrid(i,j,2)
               r = great_circle_dist( p1, p2, radius )
               if (r < r0) then
                  delp(i,j,1) = phis(i,j) + gh0*0.5*(1.0+cos(PI*r/r0))
               else
                  delp(i,j,1) = phis(i,j)
               endif
            enddo
         enddo
         initWindsCase=initWindsCase1
      case(2)
         Ubar = (2.0*pi*radius)/(12.0*86400.0)
         gh0  = 2.94e4
         phis = 0.0
         do j=js,je
            do i=is,ie
#ifdef FIVE_AVG
               pt5 = gh0 - (radius*omega*Ubar + (Ubar*Ubar)/2.) * &
                             ( -1.*cos(agrid(i  ,j  ,1))*cos(agrid(i  ,j  ,2))*sin(alpha) + &
                                   sin(agrid(i  ,j  ,2))*cos(alpha) ) ** 2.0
               pt1 = gh0 - (radius*omega*Ubar + (Ubar*Ubar)/2.) * &
                             ( -1.*cos(grid(i  ,j  ,1))*cos(grid(i  ,j  ,2))*sin(alpha) + &
                                   sin(grid(i  ,j  ,2))*cos(alpha) ) ** 2.0
               pt2 = gh0 - (radius*omega*Ubar + (Ubar*Ubar)/2.) * &
                             ( -1.*cos(grid(i+1,j  ,1))*cos(grid(i+1,j  ,2))*sin(alpha) + &
                                   sin(grid(i+1,j  ,2))*cos(alpha) ) ** 2.0
               pt3 = gh0 - (radius*omega*Ubar + (Ubar*Ubar)/2.) * &
                             ( -1.*cos(grid(i+1,j+1,1))*cos(grid(i+1,j+1,2))*sin(alpha) + &
                                   sin(grid(i+1,j+1,2))*cos(alpha) ) ** 2.0
               pt4 = gh0 - (radius*omega*Ubar + (Ubar*Ubar)/2.) * &
                             ( -1.*cos(grid(i,j+1,1))*cos(grid(i,j+1,2))*sin(alpha) + &
                                   sin(grid(i,j+1,2))*cos(alpha) ) ** 2.0
               delp(i,j,1) = (0.25*(pt1+pt2+pt3+pt4) + 3.*pt5) / 4.
#else
               delp(i,j,1) = gh0 - (radius*omega*Ubar + (Ubar*Ubar)/2.) * &
                             ( -1.*cos(agrid(i  ,j  ,1))*cos(agrid(i  ,j  ,2))*sin(alpha) + &
                                   sin(agrid(i  ,j  ,2))*cos(alpha) ) ** 2.0
#endif
            enddo
         enddo
         initWindsCase=initWindsCase2
      case(3)
!----------------------------
! Non-rotating potential flow
!----------------------------
#ifdef NO_WIND
         ubar = 0.
#else
         ubar = 40.
#endif
         gh0  = 1.0e3 * grav
         phis = 0.0
         r0 = radius/3. !RADIUS radius/3.
         p1(1) = pi*1.5
         p1(2) = 0.
         do j=jsd,jed
            do i=isd,ied
               p2(1) = agrid(i,j,1)
               p2(2) = agrid(i,j,2)
               r = great_circle_dist( p1, p2, radius )
               if (r < r0) then
                  delp(i,j,1) = phis(i,j) + gh0*0.5*(1.0+cos(PI*r/r0))
               else
                  delp(i,j,1) = phis(i,j)
               endif
! Add a constant:
               delp(i,j,1) = delp(i,j,1) + grav*2.e3
            enddo
         enddo

#ifdef NO_WIND
         u  = 0.;   v = 0.
         f0 = 0.;  fC = 0.
#else

         do j=js,je
            do i=is,ie+1
               p1(:) = grid(i  ,j ,1:2)
               p2(:) = grid(i,j+1 ,1:2)
               call mid_pt_sphere(p1, p2, p3)
               call get_unit_vector(p1, p2, e2)
               call get_latlon_vector(p3, ex, ey)
               utmp = ubar * cos(p3(2))
               vtmp = 0.
               v(i,j,1) = utmp*inner_prod(e2,ex) + vtmp*inner_prod(e2,ey)
            enddo
         enddo
         do j=js,je+1
            do i=is,ie
               p1(:) = grid(i,  j,1:2)
               p2(:) = grid(i+1,j,1:2)
               call mid_pt_sphere(p1, p2, p3)
               call get_unit_vector(p1, p2, e1)
               call get_latlon_vector(p3, ex, ey)
               utmp = ubar * cos(p3(2))
               vtmp = 0.
               u(i,j,1) = utmp*inner_prod(e1,ex) + vtmp*inner_prod(e1,ey)
            enddo
         enddo

         anti_rot = -ubar/ radius
         do j=jsd,jed+1
            do i=isd,ied+1
               fC(i,j) = 2.*anti_rot*sin(grid(i,j,2))
            enddo
         enddo
         do j=jsd,jed
            do i=isd,ied
               f0(i,j) = 2.*anti_rot*sin(agrid(i,j,2))
            enddo
         enddo
#endif
         initWindsCase= -1

      case(4)

!----------------------------
! Tropical cyclones
!----------------------------
!        f0 = 0.;  fC = 0.          ! non-rotating planet setup
          u = 0.
          v = 0.
         phis = 0.0                 ! flat terrain

         ubar = 50.                 ! maxmium wind speed (m/s)
           r0 = 250.e3              ! RADIUS of the maximum wind of the Rankine vortex
          gh0 = grav * 1.e3
 
        do j=jsd,jed
           do i=isd,ied
              delp(i,j,1) = gh0
           enddo
        enddo

!       ddeg = 2.*r0/radius     ! no merger
        ddeg = 1.80*r0/radius   ! merged 

        p1(1) = pi*1.5 - ddeg
        p1(2) = pi/18.              ! 10 N
        call rankine_vortex(ubar, r0, p1, u, v)

        p2(1) = pi*1.5 + ddeg
        p2(2) = pi/18.              ! 10 N
        call rankine_vortex(ubar, r0, p2, u, v)

#ifndef SINGULAR_VORTEX
!-----------
! Anti-pole:
!-----------
        ubar = -ubar
        call latlon2xyz(p1, e1)
        do i=1,3
           e1(i) = -e1(i)
        enddo
        call cart_to_latlon(1, e1, p3(1), p3(2))
        call rankine_vortex(ubar, r0, p3, u, v)

        call latlon2xyz(p2, e1)
        do i=1,3
           e1(i) = -e1(i)
        enddo
        call cart_to_latlon(1, e1, p4(1), p4(2))
        call rankine_vortex(ubar, r0, p4, u, v)
#endif
        call mp_update_dwinds(u, v, npx, npy, npz)
        initWindsCase=-1   ! do nothing

      case(5)

         Ubar = 20.        
         gh0  = 5960.*Grav
         phis = 0.0
         r0 = PI/9.
         p1(1) = PI/2.
         p1(2) = PI/6.
         do j=js,je
            do i=is,ie
               p2(1) = agrid(i,j,1)
               p2(2) = agrid(i,j,2)
               r = MIN(r0*r0, (p2(1)-p1(1))*(p2(1)-p1(1)) + (p2(2)-p1(2))*(p2(2)-p1(2)) )
               r = SQRT(r)
               phis(i,j) = 2000.0*Grav*(1.0-(r/r0))
            enddo
         enddo
         do j=js,je
            do i=is,ie
               delp(i,j,1) =gh0 - (radius*omega*Ubar + (Ubar*Ubar)/2.) * &
                             ( -1.*cos(agrid(i  ,j  ,1))*cos(agrid(i  ,j  ,2))*sin(alpha) + &
                                   sin(agrid(i  ,j  ,2))*cos(alpha) ) ** 2  - phis(i,j)
            enddo
         enddo
         initWindsCase=initWindsCase5
      case(6)
         gh0  = 8.E3*Grav
         R    = 4.
         omg  = 7.848E-6
         rk    = 7.848E-6
         phis = 0.0
         do j=js,je
            do i=is,ie
               A = 0.5*omg*(2.*omega+omg)*(COS(agrid(i,j,2))**2) + &
                   0.25*rk*rk*(COS(agrid(i,j,2))**(r+r)) * &
                   ( (r+1)*(COS(agrid(i,j,2))**2) + (2.*r*r-r-2.) - &
                     2.*(r*r)*COS(agrid(i,j,2))**(-2.) )
               B = (2.*(omega+omg)*rk / ((r+1)*(r+2))) * (COS(agrid(i,j,2))**r) * &
                    ( (r*r+2.*r+2.) - ((r+1.)*COS(agrid(i,j,2)))**2 )
               C = 0.25*rk*rk*(COS(agrid(i,j,2))**(2.*r)) * ( &
                   (r+1) * (COS(agrid(i,j,2))**2.) - (r+2.) )
               delp(i,j,1) =gh0 + radius*radius*(A+B*COS(r*agrid(i,j,1))+C*COS(2.*r*agrid(i,j,1)))
               delp(i,j,1) = delp(i,j,1) - phis(i,j)
            enddo
         enddo
         do j=js,je
            do i=is,ie+1
               p1(:) = grid(i  ,j ,1:2)
               p2(:) = grid(i,j+1 ,1:2)
               call mid_pt_sphere(p1, p2, p3)
               call get_unit_vector(p1, p2, e2)
               call get_latlon_vector(p3, ex, ey)
               utmp = radius*omg*cos(p3(2)) +                      &
                      radius*rk*(cos(p3(2))**(R-1))*(R*sin(p3(2))**2-cos(p3(2))**2)*cos(R*p3(1)) 
               vtmp = -radius*rk*R*sin(p3(2))*sin(R*p3(1))*cos(p3(2))**(R-1)
               v(i,j,1) = utmp*inner_prod(e2,ex) + vtmp*inner_prod(e2,ey)
            enddo
         enddo
         do j=js,je+1
            do i=is,ie
               p1(:) = grid(i,  j,1:2)
               p2(:) = grid(i+1,j,1:2)
               call mid_pt_sphere(p1, p2, p3)
               call get_unit_vector(p1, p2, e1)
               call get_latlon_vector(p3, ex, ey)
               utmp = radius*omg*cos(p3(2)) +                      &
                      radius*rk*(cos(p3(2))**(R-1))*(R*sin(p3(2))**2-cos(p3(2))**2)*cos(R*p3(1)) 
               vtmp = -radius*rk*R*sin(p3(2))*sin(R*p3(1))*cos(p3(2))**(R-1)
               u(i,j,1) = utmp*inner_prod(e1,ex) + vtmp*inner_prod(e1,ey)
            enddo
         enddo
         call mp_update_dwinds(u, v, npx, npy, npz)
         call dtoa( u, v,ua,va,npx,npy,ng)
         !call mpp_update_domains( ua, va, domain, gridtype=AGRID_PARAM)
         call atoc(ua,va,uc,vc,npx,npy,ng)
         initWindsCase=initWindsCase6
      case(7)
! Barotropically unstable jet
         gh0  = 10.E3*Grav
         phis = 0.0
         r0 = radius/12.
         p2(1) = pi/2.
         p2(2) = pi/4.
         do j=js,je
            do i=is,ie
!              ftmp = gh0
! 9-point average:
!      9  4  8
!
!      5  1  3
!          
!      6  2  7
               pt1 = gh_jet(npy, agrid(i,j,2))
               call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), pa)
               pt2 = gh_jet(npy, pa(2))
               call mid_pt_sphere(grid(i+1,j,1:2), grid(i+1,j+1,1:2), pa)
               pt3 = gh_jet(npy, pa(2))
               call mid_pt_sphere(grid(i,j+1,1:2), grid(i+1,j+1,1:2), pa)
               pt4 = gh_jet(npy, pa(2))
               call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), pa)
               pt5 = gh_jet(npy, pa(2))
               pt6 = gh_jet(npy, grid(i,  j,  2))
               pt7 = gh_jet(npy, grid(i+1,j,  2))
               pt8 = gh_jet(npy, grid(i+1,j+1,2))
               pt9 = gh_jet(npy, grid(i  ,j+1,2))
               ftmp = 0.25*pt1 + 0.125*(pt2+pt3+pt4+pt5) + 0.0625*(pt6+pt7+pt8+pt9)
#ifndef NEW_PERT
               delp(i,j,1) = ftmp + 120.*grav*cos(agrid(i,j,2)) *  &
               exp( -(3.*(agrid(i,j,1)-pi))**2 ) * exp( -(15.*(agrid(i,j,2)-pi/4.))**2 )
!              phis(i,j) = ftmp
!              delp(i,j,1) = 10.E3*grav + 120.*grav*cos(agrid(i,j,2)) *  &
!              exp( -(3.*(agrid(i,j,1)-pi))**2 ) * exp( -(15.*(agrid(i,j,2)-pi/4.))**2 )
#else
! Using great circle dist:
               p1(:) = agrid(i,j,1:2)
               delp(i,j,1) = ftmp
               r = great_circle_dist(p1, p2, radius)
               if ( r < 3.*r0 ) then
                    delp(i,j,1) = delp(i,j,1) + 1000.*grav*exp(-(r/r0)**2)
               endif
#endif
            enddo
         enddo

! v-wind:
         do j=js,je
            do i=is,ie+1
               p2(:) = grid(i,j+1,1:2)
               vv1 = u_jet(p2(2))*(ee2(2,i,j+1)*cos(p2(1)) - ee2(1,i,j+1)*sin(p2(1)))
               p1(:) = grid(i,j,1:2)
               vv3 = u_jet(p1(2))*(ee2(2,i,j)*cos(p1(1)) - ee2(1,i,j)*sin(p1(1)))
! Mid-point:
               call mid_pt_sphere(p1, p2, pa)
               vv2 = u_jet(pa(2))*(ew(2,i,j,2)*cos(pa(1)) - ew(1,i,j,2)*sin(pa(1)))
! 3-point average:
               v(i,j,1) = 0.25*(vv1 + 2.*vv2 + vv3)
!              v(i,j,1) = vv2
            enddo
         enddo
! U-wind:
         do j=js,je+1
            do i=is,ie
               p1(:) = grid(i,j,1:2)
               uu1 = u_jet(p1(2))*(ee1(2,i,j)*cos(p1(1)) - ee1(1,i,j)*sin(p1(1)))
               p2(:) = grid(i+1,j,1:2)
               uu3 = u_jet(p2(2))*(ee1(2,i+1,j)*cos(p2(1)) - ee1(1,i+1,j)*sin(p2(1)))
! Mid-point:
               call mid_pt_sphere(p1, p2, pa)
               uu2 = u_jet(pa(2))*(es(2,i,j,1)*cos(pa(1)) - es(1,i,j,1)*sin(pa(1)))
! 3-point average:
               u(i,j,1) = 0.25*(uu1 + 2.*uu2 + uu3)
!              u(i,j,1) = uu2
            enddo
         enddo
         initWindsCase=initWindsCase6  ! shouldn't do anything with this

      case(8)
!----------------------------
! Non-rotating potential flow
!----------------------------
         gh0  = 5960.*Grav
         phis = 0.0
         r0 = PI/9.
         p1(1) = PI/2.
         p1(2) = PI/6.
         do j=js,je
            do i=is,ie
               p2(1) = agrid(i,j,1)
               p2(2) = agrid(i,j,2)
               r = MIN(r0*r0, (p2(1)-p1(1))*(p2(1)-p1(1)) + (p2(2)-p1(2))*(p2(2)-p1(2)) )
               r = SQRT(r)
               phis(i,j) = 2000.0*Grav*(1.0-(r/r0))
            enddo
         enddo
         do j=js,je
            do i=is,ie
               delp(i,j,1) = gh0
            enddo
         enddo
         u  = 0.;   v = 0.
         f0 = 0.;  fC = 0.
         initWindsCase= -1

      case(9)
         jm1 = jm - 1
         DDP = PI/DBLE(jm1)
         DP  = DDP
         ll_j(1) = -0.5*PI
         do j=2,jm
            ph5  = -0.5*PI + (DBLE(j-1)-0.5)*DDP
            ll_j(j) = -0.5*PI + (DBLE(j-1)*DDP)
            sine(j) = SIN(ph5)
         enddo
         cosp( 1) =  0.
         cosp(jm) =  0.
         do j=2,jm1
            cosp(j) = (sine(j+1)-sine(j)) / DP
         enddo
         do j=2,jm
            cose(j) = 0.5 * (cosp(j-1) + cosp(j))
         enddo
         cose(1) = cose(2)
         ddeg = 180./float(jm-1)
         do j=2,jm
            deg = -90. + (float(j-1)-0.5)*ddeg
            if (deg <= 0.) then
               ll_u(j) = -10.*(deg+90.)/90.
            elseif (deg <= 60.) then
               ll_u(j) = -10. +  deg
            else
               ll_u(j) = 50. - (50./30.)* (deg - 60.)
            endif
         enddo
         ll_phi(1) = 6000. * Grav
         do j=2,jm1
            ll_phi(j)=ll_phi(j-1)  - DP*sine(j) * &
                    (radius*2.*omega + ll_u(j)/cose(j))*ll_u(j)
         enddo
         phis = 0.0
         do j=js,je
            do i=is,ie
               do jj=1,jm1
                  if ( (ll_j(jj) <= agrid(i,j,2)) .and. (agrid(i,j,2) <= ll_j(jj+1)) ) then
                     delp(i,j,1)=0.5*(ll_phi(jj)+ll_phi(jj+1))
                  endif
               enddo
            enddo
         enddo

         do j=js,je
            do i=is,ie
               if (agrid(i,j,2)*todeg <= 0.0) then
                  ua(i,j,1) = -10.*(agrid(i,j,2)*todeg + 90.)/90.
               elseif (agrid(i,j,2)*todeg <= 60.0) then
                  ua(i,j,1) = -10. + agrid(i,j,2)*todeg
               else
                  ua(i,j,1) = 50. - (50./30.)* (agrid(i,j,2)*todeg - 60.)
               endif
               va(i,j,1) = 0.0
               call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p1)
               call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p2)
               call mid_pt_sphere(grid(i+1,j,1:2), grid(i+1,j+1,1:2), p3)
               call mid_pt_sphere(grid(i,j+1,1:2), grid(i+1,j+1,1:2), p4)
               if (cubed_sphere) call rotate_winds(ua(i,j,1), va(i,j,1), p1,p2,p3,p4, agrid(i,j,1:2), 2, 1)
            enddo
         enddo

         call mpp_update_domains( ua, va, domain, gridtype=AGRID_PARAM)
         call atoc(ua,va,uc,vc,npx,npy,ng)
         call mpp_update_domains( uc, vc, domain, gridtype=CGRID_NE_PARAM)
         call fill_corners(uc, vc, npx, npy, npz, VECTOR=.true., CGRID=.true.)
         call atod(ua,va, u, v,npx,npy,ng)
         call mp_update_dwinds(u, v, npx, npy, npz)
         initWindsCase=initWindsCase9


         allocate( case9_B(isd:ied,jsd:jed) )
         call get_case9_B(case9_B)
         AofT(:) = 0.0
      end select
!--------------- end s-w cases --------------------------

! Copy 3D data for Shallow Water Tests
      do z=2,npz
         delp(:,:,z) = delp(:,:,1)
      enddo

      call mpp_update_domains( delp, domain )
      call mpp_update_domains( phis, domain )
      phi0  = delp

      call init_winds(UBar, u,v,ua,va,uc,vc, initWindsCase, npx, npy, ng, ndims, nregions)
! Copy 3D data for Shallow Water Tests
      do z=2,npz
         u(:,:,z) = u(:,:,1)
         v(:,:,z) = v(:,:,1)
      enddo

      do j=js,je
         do i=is,ie
            ps(i,j) = delp(i,j,1)
         enddo
      enddo
! -------- end s-w section ----------------------------------
#else

      if (test_case==10 .or. test_case==14) then

         alpha = 0.

   ! Initialize dry atmosphere
         q(:,:,:,:) = 3.e-6
         u(:,:,:) = 0.0
         v(:,:,:) = 0.0

       if ( test_case==14 ) then
! Aqua-planet case: mean SLP=1.E5
         phis = 0.0
         call hydro_eq(npz, is, ie, js, je, ps, phis, 1.E5,      &
                       delp, ak, bk, pt, delz, ng, .false., hybrid_z)
       else
! Initialize topography
#ifdef MARS_GCM
         gh0  = 0.*Grav
#else
         gh0  = 5960.*Grav
#endif MARS_GCM
         phis = 0.0
         r0 = PI/9.
         p1(1) = PI/4.
         p1(2) = PI/6. + (7.5/180.0)*PI
         do j=js,je
            do i=is,ie
               p2(1) = agrid(i,j,1)
               p2(2) = agrid(i,j,2)
               r = MIN(r0*r0, (p2(1)-p1(1))*(p2(1)-p1(1)) + (p2(2)-p1(2))*(p2(2)-p1(2)) )
               r = SQRT(r)
               phis(i,j) = gh0*(1.0-(r/r0))
            enddo
         enddo
         call hydro_eq(npz, is, ie, js, je, ps, phis, dry_mass,  &
                       delp, ak, bk, pt, delz, ng, mountain, hybrid_z)
       endif

      else if (test_case==11) then

#ifdef CHECK_GRID
       call pmxn(agrid, npx, npy, nregions, pmin1, pmax1, i0, j0, n0)
       if ( gid==masterproc ) write(*,*) 'A grid: Min Lon=', pmin1, 'Max lon=', pmax1
       call pmxn(agrid(isd, jsd,2), npx, npy, nregions, pmin1, pmax1, i0, j0, n0)
       if ( gid==masterproc ) write(*,*) 'A grid: Min Lat=', pmin1, 'Max lat=', pmax1
       call pmxn(grid(isd:ied,jsd:jed,1), npx, npy, nregions, pmin1, pmax1, i0, j0, n0)
       if ( gid==masterproc ) write(*,*) 'B grid: Min Lon=', pmin1, 'Max lon=', pmax1
       call pmxn(grid(isd:ied,jsd:jed,2), npx, npy, nregions, pmin1, pmax1, i0, j0, n0)
       if ( gid==masterproc ) write(*,*) 'B grid: Min Lat=', pmin1, 'Max lat=', pmax1
#endif
       call surfdrv(npx, npy, grid, agrid,   &
                    area, dx, dy, dxc, dyc, phis, gid==masterproc)
       call mpp_update_domains( phis, domain )

       if ( hybrid_z ) then
            rgrav = 1./ grav
            if( npz==32 ) then
                call compute_dz_L32( npz, ztop, dz1 )
            else
!               call mpp_error(FATAL, 'You must provide a routine for hybrid_z')
                if ( gid==masterproc ) write(*,*) 'Using const DZ'
                ztop = 45.E3           ! assuming ptop = 100.
                dz1(1) = ztop / real(npz) 
                dz1(npz) = 0.5*dz1(1)
                do z=2,npz-1
                   dz1(z) = dz1(1)
                enddo
                dz1(1) = 2.*dz1(2)
            endif

            call set_hybrid_z(is, ie, js, je, ng, npz, ztop, dz1, rgrav,  &
                              phis, ze0, delz)
!           call prt_maxmin('ZE0', ze0,  is, ie, js, je, 0, npz, 1.E-3, gid==masterproc)
!           call prt_maxmin('DZ0', delz, is, ie, js, je, 0, npz, 1.   , gid==masterproc)
       endif

! Initialize dry atmosphere
       u = 0.
       v = 0.
       q(:,:,:,:) = 0.
       q(:,:,:,1) = 3.e-6

       call hydro_eq(npz, is, ie, js, je, ps, phis, dry_mass,  &
                     delp, ak, bk, pt, delz, ng, mountain, hybrid_z)

      else if ( (test_case==12) .or. (test_case==13) ) then

         q(:,:,:,:) = 3.e-6
#ifdef TEST_TRACER
          q(:,:,:,:) = 0.
          if ( ncnst==6 ) then
              do j=js,je
                 do i=is,ie
                    q(i,j,1,1:6) = 1.
                 enddo
              enddo
              do z=1,ncnst
              do j=js,je
                 do i=is,ie
                    q(i,j,npz,z) = z
                 enddo
              enddo
              enddo
          endif
#endif
    ! Initialize surface Pressure
         ps(:,:) = 1.e5
    ! Initialize detla-P
         do z=1,npz
            do j=js,je
               do i=is,ie
                  delp(i,j,z) = ak(z+1)-ak(z) + ps(i,j)*(bk(z+1)-bk(z))
               enddo
            enddo
         enddo
    ! Setup ETA auxil variable
         eta_0 = 0.252
         do z=1,npz
            eta = 0.5*( (ak(z)+ak(z+1))/1.e5 + bk(z)+bk(z+1) )
            eta_v(z) = (eta - eta_0)*PI*0.5
         enddo
    ! Initialize winds 
         Ubar = 35.0
         r0 = 1.0
         pcen(1) = PI/9.
         pcen(2) = 2.0*PI/9. 
         if (test_case == 13) then
#ifdef ALT_PERT
             u1 = 0.0
            pt0 = 3.0
#else
             u1 = 1.0
            pt0 = 0.0
#endif
             r0 = radius/10.0
         endif
         do z=1,npz
            do j=js,je
               do i=is,ie+1
                  utmp =  Ubar * COS(eta_v(z))**(3.0/2.0) * SIN(2.0*grid(i,j+1,2))**2.0
             ! Perturbation if Case==13
                  r = great_circle_dist( pcen, grid(i,j+1,1:2), radius )
                  if (-(r/r0)**2.0 > -40.0) utmp = utmp + u1*EXP(-(r/r0)**2.0) 
                  vv1 = utmp*(ee2(2,i,j+1)*cos(grid(i,j+1,1)) - ee2(1,i,j+1)*sin(grid(i,j+1,1)))

                  utmp =  Ubar * COS(eta_v(z))**(3.0/2.0) * SIN(2.0*grid(i,j,2))**2.0
             ! Perturbation if Case==13
                  r = great_circle_dist( pcen, grid(i,j,1:2), radius )
                  if (-(r/r0)**2.0 > -40.0) utmp = utmp + u1*EXP(-(r/r0)**2.0) 
                  vv3 = utmp*(ee2(2,i,j)*cos(grid(i,j,1)) - ee2(1,i,j)*sin(grid(i,j,1)))
! Mid-point:
                  p1(:) = grid(i  ,j ,1:2)
                  p2(:) = grid(i,j+1 ,1:2)
                  call mid_pt_sphere(p1, p2, pa)
                  utmp =  Ubar * COS(eta_v(z))**(3.0/2.0) * SIN(2.0*pa(2))**2.0
             ! Perturbation if Case==13
                  r = great_circle_dist( pcen, pa, radius )
                  if (-(r/r0)**2.0 > -40.0) utmp = utmp + u1*EXP(-(r/r0)**2.0) 
                  vv2 = utmp*(ew(2,i,j,2)*cos(pa(1)) - ew(1,i,j,2)*sin(pa(1)))
! 3-point average:
                  v(i,j,z) = 0.25*(vv1 + 2.*vv2 + vv3)
               enddo
            enddo
            do j=js,je+1
               do i=is,ie
                  utmp =  Ubar * COS(eta_v(z))**(3.0/2.0) * SIN(2.0*grid(i,j,2))**2.0
             ! Perturbation if Case==13
                  r = great_circle_dist( pcen, grid(i,j,1:2), radius )
                  if (-(r/r0)**2.0 > -40.0) utmp = utmp + u1*EXP(-(r/r0)**2.0)
                  uu1 = utmp*(ee1(2,i,j)*cos(grid(i,j,1)) - ee1(1,i,j)*sin(grid(i,j,1)))

                  utmp =  Ubar * COS(eta_v(z))**(3.0/2.0) * SIN(2.0*grid(i+1,j,2))**2.0
             ! Perturbation if Case==13
                  r = great_circle_dist( pcen, grid(i+1,j,1:2), radius )
                  if (-(r/r0)**2.0 > -40.0) utmp = utmp + u1*EXP(-(r/r0)**2.0)
                  uu3 = utmp*(ee1(2,i+1,j)*cos(grid(i+1,j,1)) - ee1(1,i+1,j)*sin(grid(i+1,j,1)))
! Mid-point:
                  p1(:) = grid(i  ,j  ,1:2)
                  p2(:) = grid(i+1,j  ,1:2)
                  call mid_pt_sphere(p1, p2, pa)
                  utmp =  Ubar * COS(eta_v(z))**(3.0/2.0) * SIN(2.0*pa(2))**2.0
             ! Perturbation if Case==13
                  r = great_circle_dist( pcen, pa, radius )
                  if (-(r/r0)**2.0 > -40.0) utmp = utmp + u1*EXP(-(r/r0)**2.0)
                  uu2 = utmp*(es(2,i,j,1)*cos(pa(1)) - es(1,i,j,1)*sin(pa(1)))
! 3-point average:
                  u(i,j,z) = 0.25*(uu1 + 2.*uu2 + uu3)
               enddo
            enddo
         enddo

    ! Temperature
         eta_s = 1.0 ! Surface Level
         eta_t = 0.2 ! Tropopause
         T_0 = 288.0
         delta_T = 480000.0
         lapse_rate = 0.005
         do z=1,npz
            eta = 0.5*( (ak(z)+ak(z+1))/1.e5 + bk(z)+bk(z+1) )
        !   if (gid==masterproc) print*, z, eta
            T_mean = T_0 * eta**(RDGAS*lapse_rate/Grav)
            if (eta_t > eta) T_mean = T_mean + delta_T*(eta_t - eta)**5.0

 230  format(i4.4,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14)
            press = ptop
            do zz=1,z
               press = press + delp(is,js,zz)
            enddo
            if (gid==masterproc) write(*,230) z, eta, press/100., T_mean
            do j=js,je
               do i=is,ie
! A-grid cell center: i,j
                  pt1 = T_mean + 0.75*(eta*PI*Ubar/RDGAS)*SIN(eta_v(z))*SQRT(COS(eta_v(z))) * ( &
                              ( -2.0*(SIN(agrid(i,j,2))**6.0) *(COS(agrid(i,j,2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              2.0*Ubar*COS(eta_v(z))**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(agrid(i,j,2))**3.0)*(SIN(agrid(i,j,2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
#ifndef NO_AVG13
! 9-point average: should be 2nd order accurate for a rectangular cell
!
!      9  4  8
!
!      5  1  3
!          
!      6  2  7
!
                  call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p1)
                  pt2 = T_mean + 0.75*(eta*PI*Ubar/RDGAS)*SIN(eta_v(z))*SQRT(COS(eta_v(z))) * ( &
                              ( -2.0*(SIN(p1(2))**6.0) *(COS(p1(2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              2.0*Ubar*COS(eta_v(z))**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(p1(2))**3.0)*(SIN(p1(2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
                  call mid_pt_sphere(grid(i+1,j,1:2), grid(i+1,j+1,1:2), p1)
                  pt3 = T_mean + 0.75*(eta*PI*Ubar/RDGAS)*SIN(eta_v(z))*SQRT(COS(eta_v(z))) * ( &
                              ( -2.0*(SIN(p1(2))**6.0) *(COS(p1(2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              2.0*Ubar*COS(eta_v(z))**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(p1(2))**3.0)*(SIN(p1(2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
                  call mid_pt_sphere(grid(i,j+1,1:2), grid(i+1,j+1,1:2), p1)
                  pt4 = T_mean + 0.75*(eta*PI*Ubar/RDGAS)*SIN(eta_v(z))*SQRT(COS(eta_v(z))) * ( &
                              ( -2.0*(SIN(p1(2))**6.0) *(COS(p1(2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              2.0*Ubar*COS(eta_v(z))**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(p1(2))**3.0)*(SIN(p1(2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
                  call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p1)
                  pt5 = T_mean + 0.75*(eta*PI*Ubar/RDGAS)*SIN(eta_v(z))*SQRT(COS(eta_v(z))) * ( &
                              ( -2.0*(SIN(p1(2))**6.0) *(COS(p1(2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              2.0*Ubar*COS(eta_v(z))**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(p1(2))**3.0)*(SIN(p1(2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )

                  pt6 = T_mean + 0.75*(eta*PI*Ubar/RDGAS)*SIN(eta_v(z))*SQRT(COS(eta_v(z))) * ( &
                              ( -2.0*(SIN(grid(i,j,2))**6.0) *(COS(grid(i,j,2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              2.0*Ubar*COS(eta_v(z))**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(grid(i,j,2))**3.0)*(SIN(grid(i,j,2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
                  pt7 = T_mean + 0.75*(eta*PI*Ubar/RDGAS)*SIN(eta_v(z))*SQRT(COS(eta_v(z))) * ( &
                              ( -2.0*(SIN(grid(i+1,j,2))**6.0) *(COS(grid(i+1,j,2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              2.0*Ubar*COS(eta_v(z))**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(grid(i+1,j,2))**3.0)*(SIN(grid(i+1,j,2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
                  pt8 = T_mean + 0.75*(eta*PI*Ubar/RDGAS)*SIN(eta_v(z))*SQRT(COS(eta_v(z))) * ( &
                              ( -2.0*(SIN(grid(i+1,j+1,2))**6.0) *(COS(grid(i+1,j+1,2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              2.0*Ubar*COS(eta_v(z))**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(grid(i+1,j+1,2))**3.0)*(SIN(grid(i+1,j+1,2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
                  pt9 = T_mean + 0.75*(eta*PI*Ubar/RDGAS)*SIN(eta_v(z))*SQRT(COS(eta_v(z))) * ( &
                              ( -2.0*(SIN(grid(i,j+1,2))**6.0) *(COS(grid(i,j+1,2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              2.0*Ubar*COS(eta_v(z))**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(grid(i,j+1,2))**3.0)*(SIN(grid(i,j+1,2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
                  pt(i,j,z) = 0.25*pt1 + 0.125*(pt2+pt3+pt4+pt5) + 0.0625*(pt6+pt7+pt8+pt9)
#else
                  pt(i,j,z) = pt1
#endif

#ifdef ALT_PERT
                  r = great_circle_dist( pcen, agrid(i,j,1:2), radius )
                  if ( (r/r0)**2 < 40. ) then
                        pt(i,j,z) = pt(i,j,z) + pt0*exp(-(r/r0)**2)
                  endif
#endif
               enddo
            enddo
         enddo
         if (gid==masterproc) print*,' '
      ! Surface Geopotential
         phis(:,:)=1.e25
         do j=js,je
            do i=is,ie
               pt1 = Ubar* (COS( (eta_s-eta_0)*PI/2.0 ))**(3.0/2.0) * ( &
                              ( -2.0*(SIN(agrid(i,j,2))**6.0) *(COS(agrid(i,j,2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              Ubar*COS( (eta_s-eta_0)*PI/2.0 )**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(agrid(i,j,2))**3.0)*(SIN(agrid(i,j,2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
#ifndef NO_AVG13
! 9-point average:
!
!      9  4  8
!
!      5  1  3
!          
!      6  2  7
!
               call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p1)
               pt2 = Ubar* (COS( (eta_s-eta_0)*PI/2.0 ))**(3.0/2.0) * ( &
                           ( -2.0*(SIN(p1(2))**6.0) *(COS(p1(2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                             Ubar*COS( (eta_s-eta_0)*PI/2.0 )**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(p1(2))**3.0)*(SIN(p1(2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
               call mid_pt_sphere(grid(i+1,j,1:2), grid(i+1,j+1,1:2), p1)
               pt3 = Ubar* (COS( (eta_s-eta_0)*PI/2.0 ))**(3.0/2.0) * ( &
                           ( -2.0*(SIN(p1(2))**6.0) *(COS(p1(2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                             Ubar*COS( (eta_s-eta_0)*PI/2.0 )**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(p1(2))**3.0)*(SIN(p1(2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
               call mid_pt_sphere(grid(i,j+1,1:2), grid(i+1,j+1,1:2), p1)
               pt4 = Ubar* (COS( (eta_s-eta_0)*PI/2.0 ))**(3.0/2.0) * ( &
                           ( -2.0*(SIN(p1(2))**6.0) *(COS(p1(2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                             Ubar*COS( (eta_s-eta_0)*PI/2.0 )**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(p1(2))**3.0)*(SIN(p1(2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
               call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p1)
               pt5 = Ubar* (COS( (eta_s-eta_0)*PI/2.0 ))**(3.0/2.0) * ( &
                           ( -2.0*(SIN(p1(2))**6.0) *(COS(p1(2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                             Ubar*COS( (eta_s-eta_0)*PI/2.0 )**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(p1(2))**3.0)*(SIN(p1(2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )

               pt6 = Ubar* (COS( (eta_s-eta_0)*PI/2.0 ))**(3.0/2.0) * ( &
                              ( -2.0*(SIN(grid(i,j,2))**6.0) *(COS(grid(i,j,2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              Ubar*COS( (eta_s-eta_0)*PI/2.0 )**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(grid(i,j,2))**3.0)*(SIN(grid(i,j,2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
               pt7 = Ubar* (COS( (eta_s-eta_0)*PI/2.0 ))**(3.0/2.0) * ( &
                              ( -2.0*(SIN(grid(i+1,j,2))**6.0) *(COS(grid(i+1,j,2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              Ubar*COS( (eta_s-eta_0)*PI/2.0 )**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(grid(i+1,j,2))**3.0)*(SIN(grid(i+1,j,2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
               pt8 = Ubar* (COS( (eta_s-eta_0)*PI/2.0 ))**(3.0/2.0) * ( &
                              ( -2.0*(SIN(grid(i+1,j+1,2))**6.0) *(COS(grid(i+1,j+1,2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              Ubar*COS( (eta_s-eta_0)*PI/2.0 )**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(grid(i+1,j+1,2))**3.0)*(SIN(grid(i+1,j+1,2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
               pt9 = Ubar* (COS( (eta_s-eta_0)*PI/2.0 ))**(3.0/2.0) * ( &
                              ( -2.0*(SIN(grid(i,j+1,2))**6.0) *(COS(grid(i,j+1,2))**2.0 + 1.0/3.0) + 10.0/63.0 ) * &
                              Ubar*COS( (eta_s-eta_0)*PI/2.0 )**(3.0/2.0) + &
                              ( (8.0/5.0)*(COS(grid(i,j+1,2))**3.0)*(SIN(grid(i,j+1,2))**2.0 + 2.0/3.0) - PI/4.0 )*radius*omega )
               phis(i,j) = 0.25*pt1 + 0.125*(pt2+pt3+pt4+pt5) + 0.0625*(pt6+pt7+pt8+pt9)
#else
               phis(i,j) = pt1
#endif
            enddo
         enddo

      else if ( test_case==15 .or. test_case==19 ) then
!------------------------------------
! Non-hydrostatic 3D density current:
!------------------------------------
! C100_L64; hybrid_z = .T., make_nh = .F. ,   make_hybrid_z = .false.
! Control: npz=64;  dx = 100 m; dt = 1; n_split=10

        if ( test_case == 19 ) then
             f0(:,:) = 0.
             fC(:,:) = 0.
        endif

           phis = 0.
           u = 0.
           v = 0.
           w = 0.
          t00 = 300.
          p00 = 1.E5
          pk0 = p00**kappa
! Set up vertical coordinare with constant del-z spacing:
         ztop = 6.4E3
         ze1(    1) = ztop
         ze1(npz+1) = 0.
         do k=npz,2,-1
            ze1(k) = ze1(k+1) + ztop/real(npz)
         enddo

! Provide some room for the top layer
         ze1(1) = ztop + 1.5*ztop/real(npz)

         do j=js,je
            do i=is,ie
               ps(i,j) = p00
               pe(i,npz+1,j) = p00
               pk(i,j,npz+1) = pk0
            enddo
         enddo

         do k=npz,1,-1
            do j=js,je
               do i=is,ie
                  delz(i,j,k) = ze1(k+1) - ze1(k)
                    pk(i,j,k) = pk(i,j,k+1) + grav*delz(i,j,k)/(cp_air*t00)*pk0
                    pe(i,k,j) = pk(i,j,k)**(1./kappa)
               enddo
            enddo
         enddo

         ptop = pe(is,1,js)
         if ( gid==masterproc ) write(*,*) 'Density curent testcase: model top (mb)=', ptop/100.

         do k=1,npz+1
            do j=js,je
               do i=is,ie
                  peln(i,k,j) = log(pe(i,k,j))
                   ze0(i,j,k) = ze1(k)
               enddo
            enddo
         enddo

         do k=1,npz
            do j=js,je
               do i=is,ie
                  pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
                 delp(i,j,k) =  pe(i,k+1,j)-pe(i,k,j)
                   pt(i,j,k) = t00/pk0   ! potential temp
                enddo
            enddo
         enddo

! Perturbation: center at 3 km from the ground
         pturb = 15.
         p1(1) = pi
         p1(2) = 0.

         do k=1,npz
#ifndef STD_BUBBLE
            r0 = 0.5*(ze1(k)+ze1(k+1)) - 3.2E3
#else
            r0 = (0.5*(ze1(k)+ze1(k+1)) - 3.0E3) / 2.E3
#endif
            do j=js,je
               do i=is,ie
! Impose perturbation in potential temperature: pturb
               p2(1) = agrid(i,j,1)
               p2(2) = agrid(i,j,2)
#ifndef STD_BUBBLE
               r = great_circle_dist( p1, p2, radius )
               dist = sqrt( r**2 + r0**2 ) / 3.2E3
#else
               r = great_circle_dist( p1, p2, radius ) / 4.E3
               dist = sqrt( r**2 + r0**2 )
#endif
                  if ( dist<=1. ) then
                       q(i,j,k,1) =      pk0 * pturb/pkz(i,j,k)*(cos(pi*dist)+1.)/2.
                       pt(i,j,k) = pt(i,j,k) - pturb/pkz(i,j,k)*(cos(pi*dist)+1.)/2.
                  else
                       q(i,j,k,1) = 0.
                  endif
! Transform back to temperature:
                   pt(i,j,k) = pt(i,j,k) * pkz(i,j,k)
               enddo
            enddo
          enddo

      else if ( test_case==16 ) then

! Non-rotating:
       f0(:,:) = 0.
       fC(:,:) = 0.
! Initialize dry atmosphere
       phis = 0.
       u = 0.
       v = 0.
       p00 = 1000.E2
! Set up vertical coordinare with constant del-z spacing:
       ztop = 10.E3
       call gw_1d(npz, p00, ak, bk, ptop, ztop, ppt)

       do z=1,npz+1
          pe1(z) = ak(z) + bk(z)*p00
       enddo

       ze1(npz+1) = 0.
       do z=npz,2,-1
          ze1(z) = ze1(z+1) + ztop/real(npz)
       enddo
       ze1(1) = ztop

       if ( gid==masterproc ) write(*,*) 'Model top (pa)=', ptop

       do j=jsd,jed
          do i=isd,ied
             ps(i,j) = pe1(npz+1) 
          enddo
       enddo

       do z=1,npz+1
          do j=js,je
             do i=is,ie
                  pe(i,z,j) = pe1(z) 
                peln(i,z,j) = log(pe1(z)) 
                  pk(i,j,z) = exp(kappa*peln(i,z,j))
             enddo
          enddo
       enddo

! Horizontal shape function
       p1(1) = pi
       p1(2) = 0.
       r0 = radius / 3.
       do j=js,je
          do i=is,ie
             r = great_circle_dist( p1, agrid(i,j,1:2), radius )
             if ( r<r0 ) then
                  vort(i,j) = 0.5*(1.+cos(pi*r/r0))
             else
                  vort(i,j) = 0 
             endif
          enddo
       enddo

       q = 0.
       pk0 = p00**kappa
       pturb = 10./pk0
       do z=1,npz
          zmid = sin( 0.5*(ze1(z)+ze1(z+1))*pi/ztop )
          do j=js,je
             do i=is,ie
                 pkz(i,j,z) = (pk(i,j,z+1)-pk(i,j,z))/(kappa*(peln(i,z+1,j)-peln(i,z,j)))
                delp(i,j,z) =  pe(i,z+1,j)-pe(i,z,j)  
! Impose perturbation in potential temperature: pturb
                  pt(i,j,z) = ( ppt(z) + pturb*vort(i,j)*zmid ) * pkz(i,j,z)
                  q(i,j,z,1) = q(i,j,z,1) + vort(i,j)*zmid
             enddo
          enddo
       enddo

      elseif ( test_case==17 ) then
! Initialize dry atmosphere
       phis = 0.
       u = 0.
       v = 0.
       p00 = 1000.E2
! Set up vertical coordinare with constant del-z spacing:
       ztop = 10.E3
       call gw_1d(npz, p00, ak, bk, ptop, ztop, ppt)

       do z=1,npz+1
          pe1(z) = ak(z) + bk(z)*p00
       enddo

       ze1(npz+1) = 0.
       do z=npz,2,-1
          ze1(z) = ze1(z+1) + ztop/real(npz)
       enddo
       ze1(1) = ztop

       if ( gid==masterproc ) write(*,*) 'Model top (pa)=', ptop

       do j=jsd,jed
          do i=isd,ied
             ps(i,j) = pe1(npz+1) 
          enddo
       enddo

       do z=1,npz+1
          do j=js,je
             do i=is,ie
                  pe(i,z,j) = pe1(z) 
                peln(i,z,j) = log(pe1(z)) 
                  pk(i,j,z) = exp(kappa*peln(i,z,j))
             enddo
          enddo
       enddo

! Horizontal shape function
       p1(1) = pi
       p1(2) = pi/4.
       r0 = radius / 3.
       do j=js,je
          do i=is,ie
             r = great_circle_dist( p1, agrid(i,j,1:2), radius )
             if ( r<r0 ) then
                  vort(i,j) = 0.5*(1.+cos(pi*r/r0))
             else
                  vort(i,j) = 0 
             endif
          enddo
       enddo

         pk0 = p00**kappa
       pturb = 10./pk0
       do z=1,npz
          zmid = sin( 0.5*(ze1(z)+ze1(z+1))*pi/ztop )
          do j=js,je
             do i=is,ie
                 pkz(i,j,z) = (pk(i,j,z+1)-pk(i,j,z))/(kappa*(peln(i,z+1,j)-peln(i,z,j)))
                delp(i,j,z) =  pe(i,z+1,j)-pe(i,z,j)  
! Impose perturbation in potential temperature: pturb
                  pt(i,j,z) = ( ppt(z) + pturb*vort(i,j)*zmid ) * pkz(i,j,z)
             enddo
          enddo
       enddo

      elseif ( test_case==18 ) then
         ubar = 20.
          pt0 = 288.
         n2 = grav**2 / (cp_air*pt0)

         pcen(1) = PI/2.
         pcen(2) = PI/6.

    ! Initialize surface Pressure
         do j=js,je
            do i=is,ie
               r = great_circle_dist( pcen, agrid(i,j,1:2), radius )
               phis(i,j) = grav*2.E3*exp( -(r/1500.E3)**2 )
               ps(i,j) = 930.E2 * exp( -radius*n2*ubar/(2.*grav*grav*kappa)*(ubar/radius+2.*omega)*   &
                                       (sin(agrid(i,j,2))**2-1.) - n2/(grav*grav*kappa)*phis(i,j))
            enddo
         enddo

      do z=1,npz
            do j=js,je
               do i=is,ie
                    pt(i,j,z) = pt0
                  delp(i,j,z) = ak(z+1)-ak(z) + ps(i,j)*(bk(z+1)-bk(z))
               enddo
            enddo
! v-wind:
         do j=js,je
            do i=is,ie+1
               p1(:) = grid(i  ,j ,1:2)
               p2(:) = grid(i,j+1 ,1:2)
               call mid_pt_sphere(p1, p2, p3)
               call get_unit_vector(p1, p2, e2)
               call get_latlon_vector(p3, ex, ey)
               utmp = ubar * cos(p3(2))
               vtmp = 0.
               v(i,j,z) = utmp*inner_prod(e2,ex) + vtmp*inner_prod(e2,ey)
            enddo
         enddo

! u-wind
         do j=js,je+1
            do i=is,ie
               p1(:) = grid(i,  j,1:2)
               p2(:) = grid(i+1,j,1:2)
               call mid_pt_sphere(p1, p2, p3)
               call get_unit_vector(p1, p2, e1)
               call get_latlon_vector(p3, ex, ey)
               utmp = ubar * cos(p3(2))
               vtmp = 0.
               u(i,j,z) = utmp*inner_prod(e1,ex) + vtmp*inner_prod(e1,ey)
            enddo
         enddo
      enddo

      else if ( test_case==20 .or. test_case==21 ) then
!------------------------------------
! Non-hydrostatic 3D lee vortices
!------------------------------------
        f0(:,:) = 0.
        fC(:,:) = 0.

        if ( test_case == 20 ) then
             Ubar = 4.       ! u = Ubar * cos(lat)
             ftop = 2.0E3 * grav
        else
             Ubar = 8.       ! u = Ubar * cos(lat)
             ftop = 4.0E3 * grav
        endif

        w = 0.

         do j=js,je
            do i=is,ie+1
               p1(:) = grid(i  ,j ,1:2)
               p2(:) = grid(i,j+1 ,1:2)
               call mid_pt_sphere(p1, p2, p3)
               call get_unit_vector(p1, p2, e2)
               call get_latlon_vector(p3, ex, ey)
               utmp = ubar * cos(p3(2))
               vtmp = 0.
               v(i,j,1) = utmp*inner_prod(e2,ex) + vtmp*inner_prod(e2,ey)
            enddo
         enddo
         do j=js,je+1
            do i=is,ie
               p1(:) = grid(i,  j,1:2)
               p2(:) = grid(i+1,j,1:2)
               call mid_pt_sphere(p1, p2, p3)
               call get_unit_vector(p1, p2, e1)
               call get_latlon_vector(p3, ex, ey)
               utmp = ubar * cos(p3(2))
               vtmp = 0.
               u(i,j,1) = utmp*inner_prod(e1,ex) + vtmp*inner_prod(e1,ey)
            enddo
         enddo

! copy vertically; no wind shear
        do k=2,npz
           do j=js,je+1
              do i=is,ie
                 u(i,j,k) = u(i,j,1)
              enddo
           enddo
           do j=js,je
              do i=is,ie+1
                 v(i,j,k) = v(i,j,1)
              enddo
           enddo
        enddo

! Center of the mountain:
        p1(1) = (0.5-0.125) * pi
        p1(2) = 0.
        call latlon2xyz(p1, e1)
         uu1 =  5.0E3
         uu2 = 10.0E3
         do j=js,je
            do i=is,ie
              p2(:) = agrid(i,j,1:2)
                  r = great_circle_dist( p1, p2, radius ) 
              if ( r < pi*radius ) then
#ifdef T_ANGLE
                   call latlon2xyz(p2, e2)
! eastward vector parallel to equator
                   p3(1) = p1(1) + 0.01*pi  ! arbitrary positive number
                   p3(2) = p1(2)
                   call latlon2xyz(p3, e3)
                   e2(:) = e2(:) - e1(:)
                   e3(:) = e3(:) - e1(:)
                   call normalize_vect( e2 )
                   call normalize_vect( e3 )
! Compute angle: 0 <= acos() <= pi
                   zeta = acos( e2(1)*e3(1) + e2(2)*e3(2) + e2(3)*e3(3) )
                   if ( p2(2) <= p1(2) ) then
                        zeta = 2.*pi - zeta
                   endif
#else
                   p4(:) = p2(:) - p1(:)
                   if ( abs(p4(1)) > 1.E-12 ) then
                        zeta = asin ( p4(2) / sqrt(p4(1)**2 + p4(2)**2) ) 
                   else
                        zeta = pi/2.
                   endif
                   if ( p4(1) <= 0. ) zeta = pi - zeta
#endif
                    zeta = zeta + pi/6.
                     v1 = r/uu1 * cos( zeta )
                     v2 = r/uu2 * sin( zeta )
                   phis(i,j) = ftop / ( 1. + v1**2 + v2**2 )  
              else
                   phis(i,j) = 0.
              endif
            enddo
         enddo

       if ( hybrid_z ) then
            rgrav = 1./ grav
            if( npz==32 ) then
                call compute_dz_L32( npz, ztop, dz1 )
            else
                if ( gid==masterproc ) write(*,*) 'Using const DZ'
                ztop = 15.E3
                dz1(1) = ztop / real(npz) 
                do k=2,npz
                   dz1(k) = dz1(1)
                enddo
! Make top layer thicker
                dz1(1) = max( 1.0E3, 3.*dz1(2) )   ! min 1 km
            endif

! Re-compute ztop
             ze1(npz+1) = 0.
             do k=npz,1,-1
                ze1(k) = ze1(k+1) + dz1(k)
             enddo
             ztop = ze1(1)

            call set_hybrid_z( is, ie, js, je, ng, npz, ztop, dz1, rgrav,  &
                               phis, ze0, delz )
       else
            call mpp_error(FATAL, 'This test case is only currently setup for hybrid_z')
       endif

       do k=1,npz
          do j=js,je
             do i=is,ie
                delz(i,j,k) = ze0(i,j,k+1) - ze0(i,j,k)
             enddo
          enddo
       enddo

       p00 = 1.E5        ! mean SLP
       pk0 = p00**kappa
       t00 = 300.
       pt0 = t00/pk0
        n2 = 1.E-4
        s0 = grav*grav / (cp_air*n2) 

! For constant N2, Given z --> p
       do k=1,npz+1
          pe1(k) = p00*( (1.-s0/t00) + s0/t00*exp(-n2*ze1(k)/grav) )**(1./kappa)
       enddo

       ptop = pe1(1) 
       if ( gid==masterproc ) write(*,*) 'Lee vortex testcase: model top (mb)=', ptop/100.

! Set up fake "sigma" coordinate 
       ak(1) = pe1(1)
       bk(1) = 0.
       do k=2,npz
          bk(k) = (pe1(k) - pe1(1)) / (pe1(npz+1)-pe1(1))  ! bk == sigma
          ak(k) =  pe1(1)*(1.-bk(k)) 
       enddo                                                
       ak(npz+1) = 0.
       bk(npz+1) = 1.

! Assuming constant N
       do k=2,npz+1
          do j=js,je
             do i=is,ie
                pk(i,j,k) = pk0 - (1.-exp(-n2/grav*ze0(i,j,k))) * (grav*grav)/(n2*cp_air*pt0)
                pe(i,k,j) = pk(i,j,k) ** (1./kappa)
                peln(i,k,j) = log(pe(i,k,j)) 
             enddo
          enddo
       enddo

       do j=js,je
          do i=is,ie
               pe(i,1,j) = ptop
             peln(i,1,j) = log(pe(i,1,j)) 
               pk(i,j,1) = pe(i,1,j) ** kappa
                 ps(i,j) = pe(i,npz+1,j)
          enddo
       enddo

       do k=1,npz
          do j=js,je
             do i=is,ie
                pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
               delp(i,j,k) =  pe(i,k+1,j)-pe(i,k,j)  
                 pt(i,j,k) =  pkz(i,j,k)*grav*delz(i,j,k) / ( cp_air*(pk(i,j,k)-pk(i,j,k+1)) )
              enddo
          enddo
      enddo

      endif !test_case

      call mpp_update_domains( phis, domain )

     ftop = g_sum(phis(is:ie,js:je), is, ie, js, je, ng, area, 1)
     if(gid==masterproc) write(6,*) 'mean terrain height (m)=', ftop/grav

! The flow is initially hydrostatic
     call p_var(npz, is, ie, js, je, ptop, ptop_min, delp, delz, pt, ps,   &
                pe, peln, pk, pkz, kappa, q, ng, ncnst, dry_mass, .false., mountain, &
                moist_phys, .true., k_top, nwat)

#ifdef COLUMN_TRACER
      if( ncnst>1 ) q(:,:,:,2:ncnst) = 0.0
   ! Initialize a dummy Column Tracer
         pcen(1) = PI/9.
         pcen(2) = 2.0*PI/9.
         r0 = radius/10.0
         do z=1,npz
            do j=js,je
               do i=is,ie
                  p1(:) = grid(i  ,j ,1:2)
                  p2(:) = grid(i,j+1 ,1:2)
                  call mid_pt_sphere(p1, p2, pa)
                  call get_unit_vector(p1, p2, e2)
                  call get_latlon_vector(pa, ex, ey)
             ! Perturbation Location Case==13
                  r = great_circle_dist( pcen, pa, radius )
                  if (-(r/r0)**2.0 > -40.0) q(i,j,z,1) = EXP(-(r/r0)**2.0)
               enddo
            enddo
         enddo
#endif

#endif
    call mp_update_dwinds(u, v, npx, npy, npz)
  end subroutine init_case


  subroutine rankine_vortex(ubar, r0, p1, u, v)
!----------------------------
! Rankine vortex
!----------------------------
  real, intent(in):: ubar ! max wind (m/s)
  real, intent(in):: r0   ! Radius of max wind (m)
  real, intent(in):: p1(2)   ! center position (longitude, latitude) in radian
  real, intent(inout):: u(isd:ied,  jsd:jed+1)
  real, intent(inout):: v(isd:ied+1,jsd:jed)
! local:
  real:: p2(2), p3(2), p4(2)
  real:: e1(3), e2(3), ex(3), ey(3)
  real:: vr, r, d2, cos_p, x1, y1
  real:: utmp, vtmp
  integer i, j

! Compute u-wind
  do j=js,je+1
     do i=is,ie
        call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p2)
! shift:
        p2(1) = p2(1) - p1(1)
        cos_p = sin(p2(2))*sin(p1(2)) + cos(p2(2))*cos(p1(2))*cos(p2(1))  
        r = radius*acos(cos_p)   ! great circle distance
!       if( r<0.) call mpp_error(FATAL, 'radius negative!')
        if( r<r0 ) then
            vr = ubar*r/r0
        else
            vr = ubar*r0/r
        endif
        x1 = cos(p2(2))*sin(p2(1))
        y1 = sin(p2(2))*cos(p1(2)) - cos(p2(2))*sin(p1(2))*cos(p2(1))
        d2 = max(1.e-25, sqrt(x1**2 + y1**2))
        utmp = -vr*y1/d2
        vtmp =  vr*x1/d2
        p3(1) = grid(i,j,  1) - p1(1)
        p3(2) = grid(i,j,  2)
        p4(1) = grid(i+1,j,1) - p1(1)
        p4(2) = grid(i+1,j,2)
        call get_unit_vector(p3, p4, e1)
        call get_latlon_vector(p2, ex, ey)  ! note: p2 shifted
        u(i,j) = u(i,j) + utmp*inner_prod(e1,ex) + vtmp*inner_prod(e1,ey)
      enddo
  enddo

! Compute v-wind
  do j=js,je
     do i=is,ie+1
        call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p2)
! shift:
        p2(1) = p2(1) - p1(1)
        cos_p = sin(p2(2))*sin(p1(2)) + cos(p2(2))*cos(p1(2))*cos(p2(1))  
        r = radius*acos(cos_p)   ! great circle distance
        if( r<r0 ) then
            vr = ubar*r/r0
        else
            vr = ubar*r0/r
        endif
        x1 = cos(p2(2))*sin(p2(1))
        y1 = sin(p2(2))*cos(p1(2)) - cos(p2(2))*sin(p1(2))*cos(p2(1))
        d2 = max(1.e-25, sqrt(x1**2 + y1**2))
        utmp = -vr*y1/d2
        vtmp =  vr*x1/d2
        p3(1) = grid(i,j,  1) - p1(1)
        p3(2) = grid(i,j,  2)
        p4(1) = grid(i,j+1,1) - p1(1)
        p4(2) = grid(i,j+1,2)
        call get_unit_vector(p3, p4, e2)
        call get_latlon_vector(p2, ex, ey)  ! note: p2 shifted
        v(i,j) = v(i,j) + utmp*inner_prod(e2,ex) + vtmp*inner_prod(e2,ey)
      enddo
  enddo
  end subroutine rankine_vortex



     real function gh_jet(npy, lat_in)
     integer, intent(in):: npy
     real, intent(in):: lat_in
     real lat, lon, dp, uu
     real h0, ft
     integer j,jm

      jm = 4 * npy 
!     h0 = 10.E3
      h0 = 10.157946867E3
      dp = pi / real(jm-1)

     if ( .not. allocated(gh_table) ) then
          allocate (   gh_table(jm) )
           allocate ( lats_table(jm) )
! SP:
        gh_table(1) = grav*h0 
        lats_table(1) = -pi/2.
! Using only the mid-point for integration
      do j=2,jm
         lat = -pi/2. + (real(j-1)-0.5)*dp
         uu = u_jet(lat)
         ft = 2.*omega*sin(lat)
         gh_table(j) = gh_table(j-1) - uu*(radius*ft + tan(lat)*uu) * dp
         lats_table(j) = -pi/2. + real(j-1)*dp
      enddo
     endif

     if ( lat_in <= lats_table(1) ) then
          gh_jet = gh_table(1)
          return
     endif
     if ( lat_in >= lats_table(jm) ) then
          gh_jet = gh_table(jm)
          return
     endif

! Search:
     do j=1,jm-1
        if ( lat_in >=lats_table(j) .and. lat_in<=lats_table(j+1) ) then
             gh_jet = gh_table(j) + (gh_table(j+1)-gh_table(j))/dp * (lat_in-lats_table(j))
             return
        endif
     enddo
     end function gh_jet

     real function u_jet(lat)
      real lat, lon, dp
      real umax, en, ph0, ph1

      umax = 80.
      ph0 = pi/7.
      ph1 = pi/2. - ph0
      en =  exp( -4./(ph1-ph0)**2 )

      if ( lat>ph0 .and. lat<ph1 ) then
           u_jet = (umax/en)*exp( 1./( (lat-ph0)*(lat-ph1) ) )
      else
           u_jet = 0.
      endif
     end function u_jet
     
      subroutine get_case9_B(B)
      real, intent(OUT) :: B(isd:ied,jsd:jed)
      real :: myC,yy,myB
      integer :: i,j
! Generate B forcing function
!
      gh0 = 720.*grav
      do j=jsd,jed
         do i=isd,ied
            if (sin(agrid(i,j,2)) > 0.) then
               myC = sin(agrid(i,j,1))
                yy = (cos(agrid(i,j,2))/sin(agrid(i,j,2)))**2
               myB = gh0*yy*exp(1.-yy)
               B(i,j) = myB*myC
            else
               B(i,j) = 0.
            endif
         enddo
      enddo

   end subroutine get_case9_B
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!     
   subroutine case9_forcing1(phis,time_since_start)

   real , intent(INOUT) :: phis(isd:ied  ,jsd:jed  )
   real , intent(IN) :: time_since_start
   real :: tday, amean
   integer :: i,j
!
! Generate B forcing function
!
              tday = time_since_start/86400.0
              if (tday >= 20.) then
                 AofT(2) = 0.5*(1.-cos(0.25*PI*(tday-20)))
                 if (tday == 24) AofT(2) = 1.0
              elseif (tday <= 4.) then
                 AofT(2) = 0.5*(1.-cos(0.25*PI*tday))
              elseif (tday <= 16.) then
                 AofT(2) = 1.
              else
                 AofT(2) = 0.5*(1.+cos(0.25*PI*(tday-16.)))
              endif
              amean = 0.5*(AofT(1)+AofT(2))
              do j=jsd,jed
                 do i=isd,ied
                    phis(i,j) = amean*case9_B(i,j)
                enddo
             enddo

   end subroutine case9_forcing1
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!     
   subroutine case9_forcing2(phis)
     real ,      intent(INOUT) :: phis(isd:ied  ,jsd:jed  )
     integer :: i,j
!
! Generate B forcing function
!
          do j=jsd,jed
             do i=isd,ied
                phis(i,j) = AofT(2)*case9_B(i,j)
             enddo
          enddo
          AofT(1) = AofT(2)

   end subroutine case9_forcing2
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

      subroutine get_latlon_vector (pp, elon, elat)
      real, intent(IN)  :: pp(2)
      real, intent(OUT) :: elon(3), elat(3)

         elon(1) = -SIN(pp(1))
         elon(2) =  COS(pp(1))
         elon(3) =  0.0
         elat(1) = -SIN(pp(2))*COS(pp(1))
         elat(2) = -SIN(pp(2))*SIN(pp(1))
#ifdef RIGHT_HAND
         elat(3) =  COS(pp(2))
#else
! Left-hand system needed to be consistent with rest of the codes
         elat(3) = -COS(pp(2))
#endif

      end subroutine get_latlon_vector

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!     
!      get_stats :: get L-1, L-2, and L-inf norms and other stats as defined
!                                                in Williamson, 1994 (p.16)
       subroutine get_stats(dt, dtout, nt, maxnt, ndays, u,v,pt,delp,q,phis, ps, &
                            uc,vc, ua,va, npx, npy, npz, ncnst, ndims, nregions,    &
                            stats_lun, consv_lun, monitorFreq)
         integer,      intent(IN) :: nt, maxnt
         real  ,    intent(IN) :: dt, dtout, ndays
         real ,      intent(INOUT) ::    u(isd:ied  ,jsd:jed+1,npz)
         real ,      intent(INOUT) ::    v(isd:ied+1,jsd:jed  ,npz)
         real ,      intent(INOUT) ::   pt(isd:ied  ,jsd:jed  ,npz)
         real ,      intent(INOUT) :: delp(isd:ied  ,jsd:jed  ,npz)
         real ,      intent(INOUT) ::    q(isd:ied  ,jsd:jed  ,npz, ncnst)
         real ,      intent(INOUT) :: phis(isd:ied  ,jsd:jed  )
         real ,      intent(INOUT) ::   ps(isd:ied  ,jsd:jed  )
         real ,      intent(INOUT) ::   uc(isd:ied+1,jsd:jed  ,npz)
         real ,      intent(INOUT) ::   vc(isd:ied  ,jsd:jed+1,npz)
         real ,      intent(INOUT) ::   ua(isd:ied  ,jsd:jed  ,npz)
         real ,      intent(INOUT) ::   va(isd:ied  ,jsd:jed  ,npz)
         integer,      intent(IN) :: npx, npy, npz, ncnst
         integer,      intent(IN) :: ndims
         integer,      intent(IN) :: nregions
         integer,      intent(IN) :: stats_lun
         integer,      intent(IN) :: consv_lun
         integer,      intent(IN) :: monitorFreq

         real   :: L1_norm
         real   :: L2_norm
         real   :: Linf_norm
         real   :: pmin, pmin1, uamin1, vamin1
         real   :: pmax, pmax1, uamax1, vamax1
         real(kind=4) :: arr_r4(5)
         real   :: tmass0, tvort0, tener0, tKE0
         real   :: tmass, tvort, tener, tKE
         real   :: temp(is:ie,js:je)
         integer :: i0, j0, k0, n0
         integer :: i, j, k, n, iq

         real :: psmo, Vtx, p, w_p
         real :: x1,y1,z1,x2,y2,z2,ang

         real   :: p1(2), p2(2), p3(2), r, r0, dist, heading

         real :: uc0(isd:ied+1,jsd:jed  ,npz)
         real :: vc0(isd:ied  ,jsd:jed+1,npz)

         real :: myDay
         integer :: myRec

         myDay = ndays*((FLOAT(nt)/FLOAT(maxnt)))

#if defined(SW_DYNAMICS)
      if (test_case==0) then
         phi0 = 0.0
         do j=js,je
            do i=is,ie
               x1 = agrid(i,j,1)
               y1 = agrid(i,j,2)
               z1 = radius
               p = p0 * cos(y1)
               Vtx = ((3.0*SQRT(2.0))/2.0) * (( 1.0/cosh(p) )**2.0) * tanh(p)
               w_p = 0.0
               if (p /= 0.0) w_p = Vtx/p
              ! delp(i,j,1) = 1.0 - tanh( (p/rgamma) * sin(x1 - w_p*(nt*dt/86400.0)) )
               phi0(i,j,1) = 1.0 - tanh( (p/rgamma) * sin(x1 - w_p*(nt*dt/86400.0)) )
            enddo
         enddo
      elseif (test_case==1) then
! Get Current Height Field "Truth"
         p1(1) = pi/2.  + pi_shift
         p1(2) = 0.
         p2(1) = 3.*pi/2.  + pi_shift
         p2(2) = 0.
         r0 = radius/3. !RADIUS 3.
         dist = 2.0*pi*radius* ((FLOAT(nt)/FLOAT(maxnt)))
         heading = 3.0*pi/2.0 - alpha !5.0*pi/2.0 - alpha
         call get_pt_on_great_circle( p1, p2, dist, heading, p3)
         phi0 = 0.0
         do j=js,je
            do i=is,ie
               p2(1) = agrid(i,j,1)
               p2(2) = agrid(i,j,2)
               r = great_circle_dist( p3, p2, radius )
               if (r < r0) then
                  phi0(i,j,1) = phis(i,j) + gh0*0.5*(1.0+cos(PI*r/r0))
               else
                  phi0(i,j,1) = phis(i,j)
               endif
            enddo
         enddo
     endif

! Get Height Field Stats
         call pmxn(delp(:,:,1), npx, npy, nregions, pmin1, pmax1, i0, j0, n0)
         pmin1=pmin1/Grav
         pmax1=pmax1/Grav
         if (test_case <= 2) then
            call get_scalar_stats( delp(:,:,1), phi0(:,:,1), npx, npy, ndims, nregions, &
                                   pmin, pmax, L1_norm, L2_norm, Linf_norm)
            pmin=pmin/Grav
            pmax=pmax/Grav
            arr_r4(1) = pmin1
            arr_r4(2) = pmax1
            arr_r4(3) = L1_norm
            arr_r4(4) = L2_norm
            arr_r4(5) = Linf_norm
            if (gid == masterproc) write(stats_lun,rec=(nt)*2 + 1) arr_r4
         else
            arr_r4(1) = pmin1
            arr_r4(2) = pmax1
            arr_r4(3:5) = 0.
            pmin      = 0.
            pmax      = 0.
            L1_norm   = 0.
            L2_norm   = 0.
            Linf_norm = 0.
         endif

 200  format(i6.6,A,i6.6,A,e21.14)
 201  format('          ',A,e21.14,' ',e21.14)
 202  format('          ',A,i4.4,'x',i4.4,'x',i4.4)

         if ( (gid == masterproc) .and. MOD(nt,monitorFreq)==0 ) then
             write(*,200) nt, ' step of ', maxnt, ' DAY ', myDay
             write(*,201) 'Height MAX        : ', pmax1
             write(*,201) 'Height MIN        : ', pmin1
             write(*,202) 'HGT MAX location  : ', i0, j0, n0
             if (test_case <= 2) then
                write(*,201) 'Height L1_norm    : ', L1_norm
                write(*,201) 'Height L2_norm    : ', L2_norm
                write(*,201) 'Height Linf_norm  : ', Linf_norm
             endif
         endif

! Get UV Stats
         call dtoa(u , v , ua, va, npx, npy, ng)
         call pmxn(ua(:,:,1), npx, npy, nregions, pmin1, pmax1, i0, j0, n0)
         if (test_case <= 2) then
            call get_vector_stats( ua(:,:,1), ua0(:,:,1), va(:,:,1), va0(:,:,1), npx, npy, ndims, nregions, &
                                   pmin, pmax, L1_norm, L2_norm, Linf_norm)
         endif
         arr_r4(1) = pmin1
         arr_r4(2) = pmax1
         arr_r4(3) = L1_norm
         arr_r4(4) = L2_norm
         arr_r4(5) = Linf_norm
         if (gid == masterproc) write(stats_lun,rec=(nt)*2 + 2) arr_r4
         if ( (gid == masterproc) .and. MOD(nt,monitorFreq)==0) then
             write(*,201) 'UV     MAX        : ', pmax1
             write(*,201) 'UV     MIN        : ', pmin1
             write(*,202) 'UV  MAX location  : ', i0, j0, n0
             if (test_case <= 2) then
                write(*,201) 'UV     L1_norm    : ', L1_norm
                write(*,201) 'UV     L2_norm    : ', L2_norm
                write(*,201) 'UV     Linf_norm  : ', Linf_norm
             endif
         endif
#else

 200  format(i6.6,A,i6.6,A,e10.4)
 201  format('          ',A,e10.4,' ',e10.4,' ',i4.4,'x',i4.4,'x',i4.4,'x',i4.4)
 202  format('          ',A,e10.4,' ',e10.4,' ',i4.4,'x',i4.4,'x',i4.4,'x',i4.4,' ',e10.4)
 203  format('          ',A,i3.3,A,e10.4,' ',e10.4,' ',i4.4,'x',i4.4,'x',i4.4,'x',i4.4)

      if(gid==masterproc) write(*,200) nt, ' step of ', maxnt, ' DAY ', myDay

! Surface Pressure
     psmo = globalsum(ps(is:ie,js:je), npx, npy, is,ie, js,je)
     if(gid==masterproc) write(6,*) '         Total surface pressure =', 0.01*psmo
     call pmxn(ps, npx, npy, nregions, pmin, pmax, i0, j0, n0)
     if (gid == masterproc) then
        write(*,201) 'PS   MAX|MIN      : ', 0.01*pmax, 0.01*pmin, i0, j0, n0
     endif

! Get PT Stats
         pmax1 = -1.e25 
         pmin1 =  1.e25  
         i0=-999
         j0=-999
         k0=-999
         n0=-999
         do k=1,npz 
            call pmxn(pt(:,:,k), npx, npy, nregions, pmin, pmax, i0, j0, n0)
            pmin1 = min(pmin, pmin1)
            pmax1 = max(pmax, pmax1)
            if (pmax1 == pmax) k0 = k
         enddo
         if (gid == masterproc) then
             write(*,201) 'PT   MAX|MIN      : ', pmax1, pmin1, i0, j0, k0, n0
         endif

#if defined(DEBUG)
     if(gid==masterproc) write(6,*) ' '
         do k=1,npz
            pmax1 = -1.e25
            pmin1 =  1.e25
            i0=-999
            j0=-999
            k0=-999
            n0=-999
            call pmxn(pt(:,:,k), npx, npy, nregions, pmin, pmax, i0, j0, n0)
            pmin1 = min(pmin, pmin1)
            pmax1 = max(pmax, pmax1)
            if (gid == masterproc) then
                write(*,202) 'PT   MAX|MIN      : ', pmax1, pmin1, i0, j0, k, n0, 0.5*( (ak(k)+ak(k+1))/1.e5 + bk(k)+bk(k+1) )
            endif
         enddo
     if(gid==masterproc) write(6,*) ' '
#endif

! Get DELP Stats
         pmax1 = -1.e25 
         pmin1 =  1.e25 
         i0=-999
         j0=-999
         k0=-999
         n0=-999
         do k=1,npz
            call pmxn(delp(:,:,k), npx, npy, nregions, pmin, pmax, i0, j0, n0)
            pmin1 = min(pmin, pmin1)
            pmax1 = max(pmax, pmax1)
            if (pmax1 == pmax) k0 = k
         enddo
         if (gid == masterproc) then
             write(*,201) 'Delp MAX|MIN      : ', pmax1, pmin1, i0, j0, k0, n0
         endif

! Get UV Stats
         uamax1 = -1.e25
         uamin1 =  1.e25
         i0=-999
         j0=-999
         k0=-999
         n0=-999
         do k=1,npz
            call dtoa(u(isd,jsd,k), v(isd,jsd,k), ua(isd,jsd,k), va(isd,jsd,k), npx, npy, ng)
            call pmxn(ua(:,:,k), npx, npy, nregions, pmin, pmax, i0, j0, n0)
            uamin1 = min(pmin, uamin1)
            uamax1 = max(pmax, uamax1)
            if (uamax1 == pmax) k0 = k
         enddo
         if (gid == masterproc) then
             write(*,201) 'U    MAX|MIN      : ', uamax1, uamin1, i0, j0, k0, n0
         endif

         vamax1 = -1.e25
         vamin1 =  1.e25
         i0=-999
         j0=-999
         k0=-999
         n0=-999
         do k=1,npz
            call pmxn(va(:,:,k), npx, npy, nregions, pmin, pmax, i0, j0, n0)
            vamin1 = min(pmin, vamin1)
            vamax1 = max(pmax, vamax1)
            if (vamax1 == pmax) k0 = k
         enddo
         if (gid == masterproc) then
             write(*,201) 'V    MAX|MIN      : ', vamax1, vamin1, i0, j0, k0, n0
         endif

! Get Q Stats
         pmax1 = -1.e25 
         pmin1 =  1.e25 
         i0=-999
         j0=-999
         k0=-999
         n0=-999
         do k=1,npz
            call pmxn(q(isd,jsd,k,1), npx, npy, nregions, pmin, pmax, i0, j0, n0)
            pmin1 = min(pmin, pmin1)
            pmax1 = max(pmax, pmax1)
            if (pmax1 == pmax) k0 = k
         enddo
         if (gid == masterproc) then
             write(*,201) 'Q    MAX|MIN      : ', pmax1, pmin1, i0, j0, k0, n0
         endif

! Get tracer Stats
       do iq=2,ncnst
         pmax1 = -1.e25
         pmin1 =  1.e25
         i0=-999
         j0=-999
         k0=-999
         n0=-999
         do k=1,npz
            call pmxn(q(isd,jsd,k,iq), npx, npy, nregions, pmin, pmax, i0, j0, n0)
            pmin1 = min(pmin, pmin1)
            pmax1 = max(pmax, pmax1)
            if (pmax1 == pmax) k0 = k
         enddo
         if (gid == masterproc) then
             write(*,203) 'TR',iq-1,' MAX|MIN      : ', pmax1, pmin1, i0, j0, k0, n0
         endif
       enddo

#endif

      if (test_case == 12) then
! Get UV Stats
          call get_vector_stats( ua(:,:,22), ua0(:,:,22), va(:,:,22), va0(:,:,22), npx, npy, ndims, nregions, &
                                 pmin, pmax, L1_norm, L2_norm, Linf_norm)
          if (gid == masterproc) then
             write(*,201) 'UV(850) L1_norm    : ', L1_norm
             write(*,201) 'UV(850) L2_norm    : ', L2_norm
             write(*,201) 'UV(850) Linf_norm  : ', Linf_norm
          endif
      endif 

      tmass = 0.0
      tKE   = 0.0
      tener = 0.0
      tvort = 0.0
#if defined(SW_DYNAMICS)
      do k=1,1
#else
      do k=1,npz
#endif
! Get conservation Stats

! Conservation of Mass
         temp(:,:) = delp(is:ie,js:je,k)
         tmass0 = globalsum(temp, npx, npy, is,ie, js,je)
         tmass = tmass + tmass0

         call atoc(ua(isd,jsd,k),va(isd,jsd,k),uc0(isd,jsd,k),vc0(isd,jsd,k),npx,npy,ng)
! Conservation of Kinetic Energy
         do j=js,je
            do i=is,ie
                  temp(i,j) = ( uc0(i,j,k)*uc0(i,j,k) + uc0(i+1,j,k)*uc0(i+1,j,k) + &
                                vc0(i,j,k)*vc0(i,j,k) + vc0(i,j+1,k)*vc0(i,j+1,k) )
            enddo
         enddo
         tKE0 = globalsum(temp, npx, npy, is,ie, js,je)
         tKE = tKE + tKE0

! Conservation of Energy
         do j=js,je
            do i=is,ie
                  temp(i,j) = 0.5 * (delp(i,j,k)/Grav) * temp(i,j)  ! Include Previously calcullated KE 
                  temp(i,j) = temp(i,j) + &
                          Grav*((delp(i,j,k)/Grav + phis(i,j))*(delp(i,j,k)/Grav + phis(i,j))) - &
                          phis(i,j)*phis(i,j)
            enddo
         enddo
         tener0 = globalsum(temp, npx, npy, is,ie, js,je)
         tener = tener + tener0

! Conservation of Potential Enstrophy
         if (test_case>1) then
            do j=js,je
               do i=is,ie
                  temp(i,j) =  f0(i,j) + (1./area(i,j)) * ( (v(i+1,j,k)*dy(i+1,j) - v(i,j,k)*dy(i,j)) - &
                                                            (u(i,j+1,k)*dx(i,j+1) - u(i,j,k)*dx(i,j)) )
                  temp(i,j) = ( Grav*(temp(i,j)*temp(i,j))/delp(i,j,k) )
               enddo
            enddo
            tvort0 = globalsum(temp, npx, npy, is,ie, js,je)
            tvort = tvort + tvort0
         else
            tvort=1.
         endif
      enddo

         if (nt == 0) then
            tmass_orig = tmass
            tener_orig = tener
            tvort_orig = tvort
         endif 
         arr_r4(1) = (tmass-tmass_orig)/tmass_orig
         arr_r4(2) = (tener-tener_orig)/tener_orig
         arr_r4(3) = (tvort-tvort_orig)/tvort_orig
         arr_r4(4) = tKE
         if (test_case==12) arr_r4(4) = L2_norm 
#if defined(SW_DYNAMICS)
         myRec = nt+1
#else
         myRec = myDay*86400.0/dtout + 1 
#endif
         if (gid == masterproc) write(consv_lun,rec=myRec) arr_r4(1:4)
#if defined(SW_DYNAMICS)
         if ( (gid == masterproc) .and. MOD(nt,monitorFreq)==0) then
#else
         if ( (gid == masterproc) ) then 
#endif
             write(*,201) 'MASS TOTAL        : ', tmass
             write(*,201) 'NORMALIZED MASS   : ', (tmass-tmass_orig)/tmass_orig
             if (test_case >= 2) then
                write(*,201) 'Kinetic Energy KE : ', tKE
                write(*,201) 'ENERGY TOTAL      : ', tener
                write(*,201) 'NORMALIZED ENERGY : ', (tener-tener_orig)/tener_orig
                write(*,201) 'ENSTR TOTAL       : ', tvort
                write(*,201) 'NORMALIZED ENSTR  : ', (tvort-tvort_orig)/tvort_orig
             endif
             write(*,*) ' '
         endif
      end subroutine get_stats



   subroutine get_pt_on_great_circle(p1, p2, dist, heading, p3) 
!     get_pt_on_great_circle :: Get the mid-point on a great circle given:
!                                 -2 points (Lon/Lat) to define a great circle
!                                 -Great Cirle distance between 2 defining points
!                                 -Heading
!                              compute:
!                                 Arrival Point (Lon/Lat)

         real , intent(IN)  :: p1(2), p2(2)
         real , intent(IN)  :: dist
         real , intent(IN)  :: heading
         real , intent(OUT) :: p3(2)

         real  pha, dp

         pha = dist/radius

         p3(2) = ASIN( (COS(heading)*COS(p1(2))*SIN(pha)) + (SIN(p1(2))*COS(pha)) )
         dp = ATAN2( SIN(heading)*SIN(pha)*COS(p1(2)) , COS(pha) - SIN(p1(2))*SIN(p3(2)) )
         p3(1) = MOD( (p1(1)-pi)-dp+pi , 2.*pi ) !- pi Leave at 0 to 360

      end subroutine get_pt_on_great_circle
 

!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!      get_scalar_stats: get L-1, L-2, and L-inf norms and min/max stats as defined
!                                                in Williamson, 1994 (p.16)
!                     for any var

       subroutine get_scalar_stats(var, varT, npx, npy, ndims, nregions, &
                            vmin, vmax, L1_norm, L2_norm, Linf_norm)
         integer,      intent(IN) :: npx, npy
         integer,      intent(IN) :: ndims
         integer,      intent(IN) :: nregions
         real  ,    intent(IN) ::  var(isd:ied,jsd:jed)
         real  ,    intent(IN) :: varT(isd:ied,jsd:jed)
         real  ,   intent(OUT) :: vmin
         real  ,   intent(OUT) :: vmax
         real  ,   intent(OUT) :: L1_norm
         real  ,   intent(OUT) :: L2_norm
         real  ,   intent(OUT) :: Linf_norm

         real   :: vmean
         real   :: vvar
         real   :: vmin1
         real   :: vmax1
         real   :: pdiffmn
         real   :: pdiffmx

         real   :: varSUM, varSUM2, varMAX
         real   :: gsum
         real   :: vminT, vmaxT, vmeanT, vvarT
         integer :: i0, j0, n0

         varSUM = 0.
         varSUM2 = 0.
         varMAX = 0.
         L1_norm = 0.
         L2_norm = 0.
         Linf_norm = 0.
         vmean  = 0.
         vvar   = 0.
         vmax   = 0.
         vmin   = 0.
         pdiffmn= 0.
         pdiffmx= 0.
         vmeanT = 0.
         vvarT  = 0.
         vmaxT  = 0.
         vminT  = 0.

         vmean   = globalsum(var(is:ie,js:je) , npx, npy, is,ie, js,je)
         vmeanT  = globalsum(varT(is:ie,js:je), npx, npy, is,ie, js,je)
         vmean  = vmean  / (4.0*pi)
         vmeanT = vmeanT / (4.0*pi)

         call pmxn(var, npx, npy, nregions, vmin , vmax , i0, j0, n0)
         call pmxn(varT, npx, npy, nregions, vminT, vmaxT, i0, j0, n0)
         call pmxn(var-varT, npx, npy, nregions, pdiffmn, pdiffmx, i0, j0, n0)

         vmax = (vmax - vmaxT) / (vmaxT-vminT)
         vmin = (vmin - vminT) / (vmaxT-vminT)

         varSUM  = globalsum(varT(is:ie,js:je), npx, npy, is,ie, js,je)
         varSUM2 = globalsum(varT(is:ie,js:je)**2., npx, npy, is,ie, js,je)
         L1_norm = globalsum(ABS(var(is:ie,js:je)-varT(is:ie,js:je)), npx, npy, is,ie, js,je)
         L2_norm = globalsum((var(is:ie,js:je)-varT(is:ie,js:je))**2., npx, npy, is,ie, js,je)
         L1_norm = L1_norm/varSUM
         L2_norm = SQRT(L2_norm)/SQRT(varSUM2)

         call pmxn(ABS(varT), npx, npy, nregions, vmin, vmax, i0, j0, n0)
         varMAX = vmax
         call pmxn(ABS(var-varT), npx, npy, nregions, vmin, vmax, i0, j0, n0)
         Linf_norm = vmax/varMAX

      end subroutine get_scalar_stats
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!      get_vector_stats: get L-1, L-2, and L-inf norms and min/max stats as defined
!                                                in Williamson, 1994 (p.16)
!                     for any var

       subroutine get_vector_stats(varU, varUT, varV, varVT, &
                            npx, npy, ndims, nregions, &
                            vmin, vmax, L1_norm, L2_norm, Linf_norm)
         integer,      intent(IN) :: npx, npy
         integer,      intent(IN) :: ndims
         integer,      intent(IN) :: nregions
         real  ,    intent(IN) ::  varU(isd:ied,jsd:jed)
         real  ,    intent(IN) :: varUT(isd:ied,jsd:jed)
         real  ,    intent(IN) ::  varV(isd:ied,jsd:jed)
         real  ,    intent(IN) :: varVT(isd:ied,jsd:jed)
         real  ,   intent(OUT) :: vmin
         real  ,   intent(OUT) :: vmax
         real  ,   intent(OUT) :: L1_norm
         real  ,   intent(OUT) :: L2_norm
         real  ,   intent(OUT) :: Linf_norm

         real   ::  var(isd:ied,jsd:jed)
         real   :: varT(isd:ied,jsd:jed)
         real   :: vmean
         real   :: vvar
         real   :: vmin1
         real   :: vmax1
         real   :: pdiffmn
         real   :: pdiffmx

         real   :: varSUM, varSUM2, varMAX
         real   :: gsum
         real   :: vminT, vmaxT, vmeanT, vvarT
         integer :: i,j,n
         integer :: i0, j0, n0

         varSUM = 0.
         varSUM2 = 0.
         varMAX = 0.
         L1_norm = 0.
         L2_norm = 0.
         Linf_norm = 0.
         vmean  = 0.
         vvar   = 0.
         vmax   = 0.
         vmin   = 0.
         pdiffmn= 0.
         pdiffmx= 0.
         vmeanT = 0.
         vvarT  = 0.
         vmaxT  = 0.
         vminT  = 0.

         do j=js,je
            do i=is,ie
               var(i,j) = SQRT( (varU(i,j)-varUT(i,j))**2. + &
                                (varV(i,j)-varVT(i,j))**2. )
               varT(i,j) = SQRT( varUT(i,j)*varUT(i,j) + &
                                 varVT(i,j)*varVT(i,j) )
            enddo
         enddo
         varSUM  = globalsum(varT(is:ie,js:je), npx, npy, is,ie, js,je)
         L1_norm = globalsum(var(is:ie,js:je) , npx, npy, is,ie, js,je)
         L1_norm = L1_norm/varSUM

         call pmxn(varT, npx, npy, nregions, vmin, vmax, i0, j0, n0)
         varMAX = vmax
         call pmxn(var, npx, npy, nregions, vmin, vmax, i0, j0, n0)
         Linf_norm = vmax/varMAX

         do j=js,je
            do i=is,ie
               var(i,j) = ( (varU(i,j)-varUT(i,j))**2. + &
                            (varV(i,j)-varVT(i,j))**2. )
              varT(i,j) = ( varUT(i,j)*varUT(i,j) + &
                            varVT(i,j)*varVT(i,j) )
            enddo
         enddo
         varSUM  = globalsum(varT(is:ie,js:je), npx, npy, is,ie, js,je)
         L2_norm = globalsum(var(is:ie,js:je) , npx, npy, is,ie, js,je)
         L2_norm = SQRT(L2_norm)/SQRT(varSUM)

      end subroutine get_vector_stats
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     check_courant_numbers :: 
!
       subroutine check_courant_numbers(uc,vc, ndt, n_split, npx, npy, npz, noPrint)

       real, intent(IN) :: ndt
       integer, intent(IN) :: n_split
       integer, intent(IN) :: npx, npy, npz
       logical, OPTIONAL, intent(IN) :: noPrint
       real ,      intent(IN) ::   uc(isd:ied+1,jsd:jed  ,npz)
       real ,      intent(IN) ::   vc(isd:ied  ,jsd:jed+1,npz)
 
       real :: ideal_c=0.06
       real :: tolerance= 1.e-3
       real :: dt_inc, dt_orig 
       real   :: meanCy, minCy, maxCy, meanCx, minCx, maxCx

       real :: counter
       logical :: ideal 

       integer :: i,j,k
       real :: dt

       dt = ndt/real(n_split)

 300  format(i4.4,' ',i4.4,' ',i4.4,' ',i4.4,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14)

       dt_orig = dt
       dt_inc = 1
       ideal = .false.

       do while(.not. ideal)
       
         counter = 0
         minCy = missing
         maxCy = -1.*missing
         minCx = missing
         maxCx = -1.*missing
         meanCx = 0
         meanCy = 0
         do k=1,npz
         do j=js,je
            do i=is,ie+1
               minCx = MIN(minCx, ABS( (dt/dxc(i,j))*uc(i,j,k) ))
               maxCx = MAX(maxCx, ABS( (dt/dxc(i,j))*uc(i,j,k) ))
               meanCx = meanCx + ABS( (dt/dxc(i,j))*uc(i,j,k) )

        if (ABS( (dt/dxc(i,j))*uc(i,j,k) ) > 1.0) then
           counter = counter+1
           write(*,300) i,j,k,tile, ABS( (dt/dxc(i,j))*uc(i,j,k) ), dt, dxc(i,j), uc(i,j,k), counter 
           call exit(1)
        endif

            enddo
         enddo
         do j=js,je+1
            do i=is,ie
               minCy = MIN(minCy, ABS( (dt/dyc(i,j))*vc(i,j,k) ))
               maxCy = MAX(maxCy, ABS( (dt/dyc(i,j))*vc(i,j,k) ))
               meanCy = meanCy + ABS( (dt/dyc(i,j))*vc(i,j,k) )

        if (ABS( (dt/dyc(i,j))*vc(i,j,k) ) > 1.0) then
           counter = counter+1
           write(*,300) i,j,k,tile, ABS( (dt/dyc(i,j))*vc(i,j,k) ), dt, dyc(i,j), vc(i,j,k), counter
           call exit(1)
        endif

            enddo
         enddo
         enddo

         call mp_reduce_max(maxCx)
         call mp_reduce_max(maxCy)
         minCx = -minCx
         minCy = -minCy
         call mp_reduce_max(minCx)
         call mp_reduce_max(minCy)
         minCx = -minCx
         minCy = -minCy
         call mp_reduce_sum(meanCx)
         call mp_reduce_sum(meanCy)
         meanCx = meanCx/(6.0*DBLE(npx)*DBLE(npy-1))
         meanCy = meanCy/(6.0*DBLE(npx-1)*DBLE(npy))

         !if ( (ABS(maxCy-ideal_c) <= tolerance) .and. (ABS(maxCx-ideal_c) <= tolerance) ) then 
            ideal = .true. 
         !elseif (maxCy-ideal_c > 0) then
         !   dt = dt - dt_inc 
         !else
         !   dt = dt + dt_inc
         !endif

      enddo

         if ( (.not. present(noPrint)) .and. (gid == masterproc) ) then
            print*, ''
            print*, '--------------------------------------------'
            print*, 'Y-dir Courant number MIN  : ', minCy
            print*, 'Y-dir Courant number MAX  : ', maxCy
            print*, ''
            print*, 'X-dir Courant number MIN  : ', minCx
            print*, 'X-dir Courant number MAX  : ', maxCx
            print*, ''
            print*, 'X-dir Courant number MEAN : ', meanCx
            print*, 'Y-dir Courant number MEAN : ', meanCy
            print*, ''
            print*, 'NDT: ', ndt
            print*, 'n_split: ', n_split
            print*, 'DT: ', dt
            print*, ''
            print*, '--------------------------------------------'
            print*, ''
         endif

      end subroutine check_courant_numbers
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     pmxn :: find max and min of field p
!
      subroutine pmxn(p, npx, npy, nregions, pmin, pmax, i0, j0, n0)
         integer,      intent(IN) :: npx
         integer,      intent(IN) :: npy
         integer,      intent(IN) :: nregions
         real  , intent(IN)  :: p(isd:ied,jsd:jed)
         real  , intent(OUT) :: pmin
         real  , intent(OUT) :: pmax
         integer,      intent(OUT) :: i0
         integer,      intent(OUT) :: j0
         integer,      intent(OUT) :: n0

         real   :: temp
         integer :: i,j,n

         pmax = -1.e25
         pmin =  1.e25 
         i0 = -999
         j0 = -999
         n0 = tile

            do j=js,je
               do i=is,ie
                  temp = p(i,j)
                  if (temp > pmax) then
                     pmax = temp
                     i0 = i
                     j0 = j
                  elseif (temp < pmin) then
                     pmin = temp
                  endif
            enddo
         enddo

         temp = pmax
         call mp_reduce_max(temp)
         if (temp /= pmax) then
            i0 = -999
            j0 = -999
            n0 = -999
         endif
         pmax = temp
         call mp_reduce_max(i0)
         call mp_reduce_max(j0)
         call mp_reduce_max(n0)

         pmin = -pmin                  
         call mp_reduce_max(pmin)
         pmin = -pmin

      end subroutine pmxn
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     output_ncdf :: write out NETCDF fields
!
      subroutine output_ncdf(dt, nt, maxnt, nout, u,v,pt,delp,q,phis,ps, uc,vc, ua,va, &
                        omga, npx, npy, npz, ng, ncnst, ndims, nregions, ncid, &
                        npx_p1_id, npy_p1_id, npx_id, npy_id, npz_id, ntiles_id, ncnst_id, nt_id, &
                        phis_id, delp_id, ps_id, pt_id, pv_id, om_id, u_id, v_id, q_id, tracers_ids,  &
                        lats_id, lons_id)
      real,         intent(IN) :: dt
      integer,      intent(IN) :: nt, maxnt
      integer,      intent(INOUT) :: nout

      real ,      intent(INOUT) ::    u(isd:ied  ,jsd:jed+1,npz)
      real ,      intent(INOUT) ::    v(isd:ied+1,jsd:jed  ,npz)
      real ,      intent(INOUT) ::   pt(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(INOUT) :: delp(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(INOUT) ::    q(isd:ied  ,jsd:jed  ,npz, ncnst)

      real ,      intent(INOUT) :: phis(isd:ied  ,jsd:jed  )
      real ,      intent(INOUT) ::   ps(isd:ied  ,jsd:jed  )

      real ,      intent(INOUT) ::   uc(isd:ied+1,jsd:jed  ,npz)
      real ,      intent(INOUT) ::   vc(isd:ied  ,jsd:jed+1,npz)
      real ,      intent(INOUT) ::   ua(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(INOUT) ::   va(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(INOUT) :: omga(isd:ied  ,jsd:jed  ,npz)

      integer,      intent(IN) :: npx, npy, npz
      integer,      intent(IN) :: ng, ncnst
      integer,      intent(IN) :: ndims
      integer,      intent(IN) :: nregions
      integer,      intent(IN) :: ncid
      integer,      intent(IN) :: npx_p1_id, npy_p1_id, npx_id, npy_id, npz_id, ncnst_id
      integer,      intent(IN) :: ntiles_id, nt_id
      integer,      intent(IN) :: phis_id, delp_id, ps_id, pt_id, pv_id, u_id, v_id, q_id
      integer,      intent(IN) :: om_id          ! omega (dp/dt)
      integer,      intent(IN) :: tracers_ids(ncnst-1)
      integer,      intent(IN) :: lats_id, lons_id

      real, allocatable :: tmp(:,:,:)
      real, allocatable :: tmpA(:,:,:)
#if defined(SW_DYNAMICS) 
      real, allocatable :: ut(:,:,:)
      real, allocatable :: vt(:,:,:)
#else       
      real, allocatable :: ut(:,:,:,:)
      real, allocatable :: vt(:,:,:,:)
      real, allocatable :: tmpA_3d(:,:,:,:)
#endif
      real, allocatable :: vort(:,:)

      real   :: p1(2)      ! Temporary Point
      real   :: p2(2)      ! Temporary Point
      real   :: p3(2)      ! Temporary Point
      real   :: p4(2)      ! Temporary Point
      real   :: pa(2)      ! Temporary Point
      real   :: utmp, vtmp, r, r0, dist, heading
      integer   ::  i,j,k,n,iq,nreg

      real :: Vtx, p, w_p
      real :: x1,y1,z1,x2,y2,z2,ang

      allocate( tmp(npx  ,npy  ,nregions) )
      allocate( tmpA(npx-1,npy-1,nregions) )
#if defined(SW_DYNAMICS) 
      allocate( ut(npx-1,npy-1,nregions) )
      allocate( vt(npx-1,npy-1,nregions) )
#else
      allocate( ut(npx-1,npy-1,npz,nregions) )
      allocate( vt(npx-1,npy-1,npz,nregions) )
      allocate( tmpA_3d(npx-1,npy-1,npz,nregions) )
#endif
      allocate( vort(isd:ied,jsd:jed) ) 

      nout = nout + 1

      if (nt==0) then
         tmp(is:ie+1,js:je+1,tile) = grid(is:ie+1,js:je+1,2)
         call wrtvar_ncdf(ncid, lats_id, nout, is,ie+1, js,je+1, npx+1, npy+1, 1, nregions, tmp(1:npx,1:npy,1:nregions), 3)
         tmp(is:ie+1,js:je+1,tile) = grid(is:ie+1,js:je+1,1)
         call wrtvar_ncdf(ncid, lons_id, nout, is,ie+1, js,je+1, npx+1, npy+1, 1, nregions, tmp(1:npx,1:npy,1:nregions), 3)
      endif

#if defined(SW_DYNAMICS)
      if (test_case > 1) then
         tmpA(is:ie,js:je,tile) = delp(is:ie,js:je,1)/Grav

         if ((nt==0) .and. (test_case==2)) then
         Ubar = (2.0*pi*radius)/(12.0*86400.0)
         gh0  = 2.94e4
         phis = 0.0
         do j=js,je+1
            do i=is,ie+1
               tmp(i,j,tile) = (gh0 - (radius*omega*Ubar + (Ubar*Ubar)/2.) * &
                           ( -1.*cos(grid(i  ,j  ,1))*cos(grid(i  ,j  ,2))*sin(alpha) + &
                                 sin(grid(i  ,j  ,2))*cos(alpha) ) ** 2.0) / Grav
            enddo
         enddo
         endif

      else

       if (test_case==1) then
! Get Current Height Field "Truth"
         p1(1) = pi/2. + pi_shift
         p1(2) = 0.
         p2(1) = 3.*pi/2. + pi_shift
         p2(2) = 0.
         r0 = radius/3. !RADIUS /3.
         dist = 2.0*pi*radius* ((FLOAT(nt)/FLOAT(maxnt)))
         heading = 5.0*pi/2.0 - alpha
         call get_pt_on_great_circle( p1, p2, dist, heading, p3)
            do j=jsd,jed
               do i=isd,ied
                  p2(1) = agrid(i,j,1)
                  p2(2) = agrid(i,j,2)
                  r = great_circle_dist( p3, p2, radius )
                  if (r < r0) then
                     phi0(i,j,1) = phis(i,j) + gh0*0.5*(1.0+cos(PI*r/r0))
                  else
                     phi0(i,j,1) = phis(i,j)
                  endif
               enddo
            enddo
         elseif (test_case == 0) then
           phi0 = 0.0
           do j=jsd,jed
              do i=isd,ied
               x1 = agrid(i,j,1)
               y1 = agrid(i,j,2)
               z1 = radius
               p = p0 * cos(y1)
               Vtx = ((3.0*SQRT(2.0))/2.0) * (( 1.0/cosh(p) )**2.0) * tanh(p)
               w_p = 0.0
               if (p /= 0.0) w_p = Vtx/p
               phi0(i,j,1) = 1.0 - tanh( (p/rgamma) * sin(x1 - w_p*(nt*dt/86400.0)) )
              enddo
           enddo
         endif

         tmpA(is:ie,js:je,tile) = phi0(is:ie,js:je,1)
         call wrtvar_ncdf(ncid, phis_id, nout, is,ie, js,je, npx, npy, npz, nregions, tmpA, 3)
         tmpA(is:ie,js:je,tile) = delp(is:ie,js:je,1)
      endif
      call wrtvar_ncdf(ncid, ps_id, nout, is,ie, js,je, npx, npy, npz, nregions, tmpA, 3)

      if (test_case == 9) then
! Calc Vorticity
         do j=jsd,jed
            do i=isd,ied
               vort(i,j) = f0(i,j) + (1./area(i,j)) * ( (v(i+1,j,1)*dy(i+1,j) - v(i,j,1)*dy(i,j)) - &
                                                        (u(i,j+1,1)*dx(i,j+1) - u(i,j,1)*dx(i,j)) )
               vort(i,j) = Grav*vort(i,j)/delp(i,j,1)
            enddo
         enddo
         tmpA(is:ie,js:je,tile) = vort(is:ie,js:je)
         call wrtvar_ncdf(ncid, pv_id, nout, is,ie, js,je, npx, npy, npz, nregions, tmpA, 3)
      endif

      call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, 1)
      do j=js,je
         do i=is,ie
            ut(i,j,tile) = ua(i,j,1)
            vt(i,j,tile) = va(i,j,1)
         enddo
      enddo

      call wrtvar_ncdf(ncid, u_id, nout, is,ie, js,je, npx, npy, npz, nregions, ut(1:npx-1,1:npy-1,1:nregions), 3)
      call wrtvar_ncdf(ncid, v_id, nout, is,ie, js,je, npx, npy, npz, nregions, vt(1:npx-1,1:npy-1,1:nregions), 3)

      if ((test_case >= 2) .and. (nt==0) ) then
         tmpA(is:ie,js:je,tile) = phis(is:ie,js:je)/Grav
         call wrtvar_ncdf(ncid, phis_id, nout, is,ie, js,je, npx, npy, npz, nregions, tmpA, 3)
      endif
#else

! Write Moisture Data
      tmpA_3d(is:ie,js:je,1:npz,tile) = q(is:ie,js:je,1:npz,1)
      call wrtvar_ncdf(ncid, q_id, nout, is,ie, js,je, npx, npy, npz, nregions, tmpA_3d, 4)

! Write Tracer Data
      do iq=2,ncnst
         tmpA_3d(is:ie,js:je,1:npz,tile) = q(is:ie,js:je,1:npz,iq)
         call wrtvar_ncdf(ncid, tracers_ids(iq-1), nout, is,ie, js,je, npx, npy, npz, nregions, tmpA_3d, 4)
      enddo

! Write Surface height data
      tmpA(is:ie,js:je,tile) = phis(is:ie,js:je)/Grav
      call wrtvar_ncdf(ncid, phis_id, nout, is,ie, js,je, npx, npy, 1, nregions, tmpA, 3)

! Write Pressure Data
      tmpA(is:ie,js:je,tile) = ps(is:ie,js:je)
      call wrtvar_ncdf(ncid, ps_id, nout, is,ie, js,je, npx, npy, 1, nregions, tmpA, 3)
      do k=1,npz
         tmpA_3d(is:ie,js:je,k,tile) = delp(is:ie,js:je,k)/Grav
      enddo
      call wrtvar_ncdf(ncid, delp_id, nout, is,ie, js,je, npx, npy, npz, nregions, tmpA_3d, 4)

! Write PT Data
      do k=1,npz
         tmpA_3d(is:ie,js:je,k,tile) = pt(is:ie,js:je,k)
      enddo
      call wrtvar_ncdf(ncid, pt_id, nout, is,ie, js,je, npx, npy, npz, nregions, tmpA_3d, 4)

! Write U,V Data
      call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, npz, 1)
      do k=1,npz
         do j=js,je
            do i=is,ie
               ut(i,j,k,tile) = ua(i,j,k)
               vt(i,j,k,tile) = va(i,j,k)
            enddo
         enddo
      enddo
      call wrtvar_ncdf(ncid, u_id, nout, is,ie, js,je, npx, npy, npz, nregions, ut(1:npx-1,1:npy-1,1:npz,1:nregions), 4)
      call wrtvar_ncdf(ncid, v_id, nout, is,ie, js,je, npx, npy, npz, nregions, vt(1:npx-1,1:npy-1,1:npz,1:nregions), 4)


! Calc Vorticity
      do k=1,npz
         do j=js,je
            do i=is,ie
               tmpA_3d(i,j,k,tile) = rarea(i,j) * ( (v(i+1,j,k)*dy(i+1,j) - v(i,j,k)*dy(i,j)) - &
                                                    (u(i,j+1,k)*dx(i,j+1) - u(i,j,k)*dx(i,j)) )
            enddo
         enddo
      enddo
      call wrtvar_ncdf(ncid, pv_id, nout, is,ie, js,je, npx, npy, npz, nregions, tmpA_3d, 4)
!
! Output omega (dp/dt):
      do k=1,npz
         do j=js,je
            do i=is,ie
               tmpA_3d(i,j,k,tile) = omga(i,j,k)
            enddo
         enddo
      enddo
      call wrtvar_ncdf(ncid, om_id, nout, is,ie, js,je, npx, npy, npz, nregions, tmpA_3d, 4)

#endif

      deallocate( tmp )
      deallocate( tmpA )
#if defined(SW_DYNAMICS) 
      deallocate( ut )
      deallocate( vt )
#else
      deallocate( ut )
      deallocate( vt )
      deallocate( tmpA_3d )
#endif
      deallocate( vort )

      end subroutine output_ncdf
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     output :: write out fields
!
      subroutine output(dt, nt, maxnt, nout, u,v,pt,delp,q,phis,ps, uc,vc, ua,va, &
                        npx, npy, npz, ng, ncnst, ndims, nregions, phis_lun, phi_lun, &
                        pt_lun, pv_lun, uv_lun)

      real,         intent(IN) :: dt
      integer,      intent(IN) :: nt, maxnt
      integer,      intent(INOUT) :: nout

      real ,      intent(INOUT) ::    u(isd:ied  ,jsd:jed+1,npz)
      real ,      intent(INOUT) ::    v(isd:ied+1,jsd:jed  ,npz)
      real ,      intent(INOUT) ::   pt(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(INOUT) :: delp(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(INOUT) ::    q(isd:ied  ,jsd:jed  ,npz, ncnst)

      real ,      intent(INOUT) :: phis(isd:ied  ,jsd:jed  )
      real ,      intent(INOUT) ::   ps(isd:ied  ,jsd:jed  )

      real ,      intent(INOUT) ::   uc(isd:ied+1,jsd:jed  ,npz)
      real ,      intent(INOUT) ::   vc(isd:ied  ,jsd:jed+1,npz)
      real ,      intent(INOUT) ::   ua(isd:ied  ,jsd:jed  ,npz)
      real ,      intent(INOUT) ::   va(isd:ied  ,jsd:jed  ,npz)

      integer,      intent(IN) :: npx, npy, npz
      integer,      intent(IN) :: ng, ncnst
      integer,      intent(IN) :: ndims
      integer,      intent(IN) :: nregions
      integer,      intent(IN) :: phis_lun, phi_lun, pt_lun, pv_lun, uv_lun

      real   ::  tmp(1-ng:npx  +ng,1-ng:npy  +ng,1:nregions)
      real   :: tmpA(1-ng:npx-1+ng,1-ng:npy-1+ng,1:nregions)
      real   :: p1(2)      ! Temporary Point
      real   :: p2(2)      ! Temporary Point
      real   :: p3(2)      ! Temporary Point
      real   :: p4(2)      ! Temporary Point
      real   :: pa(2)      ! Temporary Point
      real   :: ut(1:npx,1:npy,1:nregions)
      real   :: vt(1:npx,1:npy,1:nregions)
      real   :: utmp, vtmp, r, r0, dist, heading
      integer   ::  i,j,k,n,nreg
      real   :: vort(isd:ied,jsd:jed)

      real :: Vtx, p, w_p
      real :: x1,y1,z1,x2,y2,z2,ang

      nout = nout + 1

#if defined(SW_DYNAMICS)
      if (test_case > 1) then
         call atob_s(delp(:,:,1)/Grav, tmp(isd:ied+1,jsd:jed+1,tile), npx,npy) !, altInterp=1)
         tmpA(is:ie,js:je,tile) = delp(is:ie,js:je,1)/Grav

         if ((nt==0) .and. (test_case==2)) then
         Ubar = (2.0*pi*radius)/(12.0*86400.0)
         gh0  = 2.94e4
         phis = 0.0
         do j=js,je+1
            do i=is,ie+1
               tmp(i,j,tile) = (gh0 - (radius*omega*Ubar + (Ubar*Ubar)/2.) * &
                           ( -1.*cos(grid(i  ,j  ,1))*cos(grid(i  ,j  ,2))*sin(alpha) + &
                                 sin(grid(i  ,j  ,2))*cos(alpha) ) ** 2.0) / Grav
            enddo
         enddo
         endif

      else

       if (test_case==1) then
! Get Current Height Field "Truth"
         p1(1) = pi/2. + pi_shift
         p1(2) = 0.
         p2(1) = 3.*pi/2. + pi_shift
         p2(2) = 0.
         r0 = radius/3. !RADIUS /3.
         dist = 2.0*pi*radius* ((FLOAT(nt)/FLOAT(maxnt)))
         heading = 5.0*pi/2.0 - alpha
         call get_pt_on_great_circle( p1, p2, dist, heading, p3)
            do j=jsd,jed
               do i=isd,ied
                  p2(1) = agrid(i,j,1)
                  p2(2) = agrid(i,j,2)
                  r = great_circle_dist( p3, p2, radius )
                  if (r < r0) then
                     phi0(i,j,1) = phis(i,j) + gh0*0.5*(1.0+cos(PI*r/r0))
                  else
                     phi0(i,j,1) = phis(i,j)
                  endif
               enddo
            enddo
         elseif (test_case == 0) then
           phi0 = 0.0
           do j=jsd,jed
              do i=isd,ied
               x1 = agrid(i,j,1) 
               y1 = agrid(i,j,2)
               z1 = radius
               p = p0 * cos(y1)
               Vtx = ((3.0*SQRT(2.0))/2.0) * (( 1.0/cosh(p) )**2.0) * tanh(p)
               w_p = 0.0
               if (p /= 0.0) w_p = Vtx/p
               phi0(i,j,1) = 1.0 - tanh( (p/rgamma) * sin(x1 - w_p*(nt*dt/86400.0)) )
              enddo
           enddo
         endif

         call atob_s(phi0(:,:,1), tmp(isd:ied+1,jsd:jed+1,tile), npx,npy) !, altInterp=1)
         tmpA(is:ie,js:je,tile) = phi0(is:ie,js:je,1)
         call wrt2d(phis_lun, nout  , is,ie, js,je, npx, npy, nregions, tmpA(1:npx-1,1:npy-1,1:nregions))
         call atob_s(delp(:,:,1), tmp(isd:ied+1,jsd:jed+1,tile), npx,npy) !, altInterp=1)
         tmpA(is:ie,js:je,tile) = delp(is:ie,js:je,1)
      endif
   !   call wrt2d(phi_lun, nout, is,ie+1, js,je+1, npx+1, npy+1, nregions, tmp(1:npx,1:npy,1:nregions))
      call wrt2d(phi_lun, nout, is,ie, js,je, npx, npy, nregions, tmpA(1:npx-1,1:npy-1,1:nregions))

      if (test_case == 9) then
! Calc Vorticity
         do j=jsd,jed
            do i=isd,ied
               vort(i,j) = f0(i,j) + (1./area(i,j)) * ( (v(i+1,j,1)*dy(i+1,j) - v(i,j,1)*dy(i,j)) - &
                                                        (u(i,j+1,1)*dx(i,j+1) - u(i,j,1)*dx(i,j)) )
               vort(i,j) = Grav*vort(i,j)/delp(i,j,1)
            enddo
         enddo
         call atob_s(vort, tmp(isd:ied+1,jsd:jed+1,tile), npx,npy) !, altInterp=1)
         call wrt2d(pv_lun, nout, is,ie+1, js,je+1, npx+1, npy+1, nregions, tmp(1:npx,1:npy,1:nregions))
      endif

      call dtoa(u , v , ua, va, npx, npy, ng)
! Rotate winds to standard Lat-Lon orientation
      if (cubed_sphere) then
         do j=js,je
            do i=is,ie
               call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p1)
               call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p2)
               call mid_pt_sphere(grid(i+1,j,1:2), grid(i+1,j+1,1:2), p3)
               call mid_pt_sphere(grid(i,j+1,1:2), grid(i+1,j+1,1:2), p4)
               utmp = ua(i,j,1)
               vtmp = va(i,j,1)
               if (cubed_sphere) call rotate_winds(utmp,vtmp, p1,p2,p3,p4, agrid(i,j,1:2), 2, 2)
               ut(i,j,tile) = utmp
               vt(i,j,tile) = vtmp
            enddo
         enddo
      endif

      call wrt2d(uv_lun, 2*(nout-1) + 1, is,ie, js,je, npx, npy, nregions,   ut(1:npx-1,1:npy-1,1:nregions))
      call wrt2d(uv_lun, 2*(nout-1) + 2, is,ie, js,je, npx, npy, nregions,   vt(1:npx-1,1:npy-1,1:nregions))

      if ((test_case >= 2) .and. (nt==0) ) then
         call atob_s(phis/Grav, tmp(isd:ied+1,jsd:jed+1,tile), npx,npy) !, altInterp=1)
       !  call wrt2d(phis_lun, nout  , is,ie+1, js,je+1, npx+1, npy+1, nregions, tmp(1:npx,1:npy,1:nregions))
         tmpA(is:ie,js:je,tile) = phis(is:ie,js:je)/Grav
         call wrt2d(phis_lun, nout  , is,ie, js,je, npx, npy, nregions, tmpA(1:npx-1,1:npy-1,1:nregions))
      endif
#else

! Write Surface height data
      if (nt==0) then
         tmpA(is:ie,js:je,tile) = phis(is:ie,js:je)/Grav
         call wrt2d(phis_lun, nout  , is,ie, js,je, npx, npy, nregions, tmpA(1:npx-1,1:npy-1,1:nregions))
      endif

! Write Pressure Data

      !if (tile==2) then
      !   do i=is,ie
      !      print*, i, ps(i,35) 
      !   enddo
      !endif
      tmpA(is:ie,js:je,tile) = ps(is:ie,js:je)
      call wrt2d(phi_lun, (nout-1)*(npz+1) + 1, is,ie, js,je, npx, npy, nregions, tmpA(1:npx-1,1:npy-1,1:nregions))
      do k=1,npz
         tmpA(is:ie,js:je,tile) = delp(is:ie,js:je,k)/Grav
         call wrt2d(phi_lun, (nout-1)*(npz+1) + 1 + k, is,ie, js,je, npx, npy, nregions, tmpA(1:npx-1,1:npy-1,1:nregions))
      enddo

! Write PT Data
      do k=1,npz
         tmpA(is:ie,js:je,tile) = pt(is:ie,js:je,k)
         call wrt2d(pt_lun, (nout-1)*npz + (k-1) + 1, is,ie, js,je, npx, npy, nregions, tmpA(1:npx-1,1:npy-1,1:nregions))
      enddo

! Write U,V Data
      do k=1,npz
         call dtoa(u(isd,jsd,k), v(isd,jsd,k), ua(isd,jsd,k), va(isd,jsd,k), npx, npy, ng)
! Rotate winds to standard Lat-Lon orientation
         if (cubed_sphere) then
            do j=js,je
               do i=is,ie
                 call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p1)
                 call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p2)
                 call mid_pt_sphere(grid(i+1,j,1:2), grid(i+1,j+1,1:2), p3)
                 call mid_pt_sphere(grid(i,j+1,1:2), grid(i+1,j+1,1:2), p4)
                 utmp = ua(i,j,k)
                 vtmp = va(i,j,k)
                 if (cubed_sphere) call rotate_winds(utmp,vtmp, p1,p2,p3,p4, agrid(i,j,1:2), 2, 2)
                 ut(i,j,tile) = utmp
                 vt(i,j,tile) = vtmp
               enddo
            enddo
         endif
         call wrt2d(uv_lun, 2*((nout-1)*npz + (k-1)) + 1, is,ie, js,je, npx, npy, nregions,   ut(1:npx-1,1:npy-1,1:nregions))
         call wrt2d(uv_lun, 2*((nout-1)*npz + (k-1)) + 2, is,ie, js,je, npx, npy, nregions,   vt(1:npx-1,1:npy-1,1:nregions))
      enddo
#endif
      end subroutine output
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!     wrt2d_ncdf :: write out a 2d field
!        
      subroutine wrtvar_ncdf(ncid, varid, nrec, i1,i2, j1,j2, npx, npy, npz, ntiles, p, ndims)
#include <netcdf.inc>
         integer,      intent(IN) :: ncid, varid
         integer,      intent(IN) :: nrec
         integer,      intent(IN) :: i1,i2,j1,j2
         integer,      intent(IN) :: npx
         integer,      intent(IN) :: npy
         integer,      intent(IN) :: npz
         integer,      intent(IN) :: ntiles
         real  , intent(IN)  :: p(npx-1,npy-1,npz,ntiles)
         integer,      intent(IN) :: ndims

         integer :: error
         real(kind=4), allocatable :: p_R4(:,:,:,:)
         integer :: i,j,k,n
         integer :: istart(ndims+1), icount(ndims+1)

         allocate( p_R4(npx-1,npy-1,npz,ntiles) )

         p_R4(:,:,:,:) = missing
         p_R4(i1:i2,j1:j2,1:npz,tile) = p(i1:i2,j1:j2,1:npz,tile)
         call mp_gather(p_R4, i1,i2, j1,j2, npx-1, npy-1, npz, ntiles)

         istart(:) = 1
         istart(ndims+1) = nrec
         icount(1) = npx-1
         icount(2) = npy-1
         icount(3) = npz
         if (ndims == 3) icount(3) = ntiles
         if (ndims == 4) icount(4) = ntiles
         icount(ndims+1) = 1

         if (gid == masterproc) then  
            error = NF_PUT_VARA_REAL(ncid, varid, istart, icount, p_R4)
         endif ! masterproc

         deallocate( p_R4 )

      end subroutine wrtvar_ncdf
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!     wrt2d :: write out a 2d field
!
      subroutine wrt2d(iout, nrec, i1,i2, j1,j2, npx, npy, nregions, p)
         integer,      intent(IN) :: iout
         integer,      intent(IN) :: nrec
         integer,      intent(IN) :: i1,i2,j1,j2
         integer,      intent(IN) :: npx
         integer,      intent(IN) :: npy
         integer,      intent(IN) :: nregions
         real  , intent(IN)  :: p(npx-1,npy-1,nregions)

         real(kind=4) :: p_R4(npx-1,npy-1,nregions)
         integer :: i,j,n

         do n=tile,tile
            do j=j1,j2
               do i=i1,i2
                  p_R4(i,j,n) = p(i,j,n)
               enddo
            enddo
         enddo

         call mp_gather(p_R4, i1,i2, j1,j2, npx-1, npy-1, nregions) 

         if (gid == masterproc) then
            write(iout,rec=nrec) p_R4(1:npx-1,1:npy-1,1:nregions)
         endif ! masterproc

      end subroutine wrt2d
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!     init_double_periodic
!
      subroutine init_double_periodic(u,v,w,pt,delp,q,phis, ps,pe,peln,pk,pkz,  uc,vc, ua,va, ak, bk,  &
                                      npx, npy, npz, ng, ncnst, nwat, k_top, ndims, nregions, dry_mass, &
                                      mountain, moist_phys, hydrostatic, hybrid_z, delz, ze0)

        real ,      intent(INOUT) ::    u(isd:ied  ,jsd:jed+1,npz)
        real ,      intent(INOUT) ::    v(isd:ied+1,jsd:jed  ,npz)
        real ,      intent(INOUT) ::    w(isd:ied  ,jsd:jed  ,npz)
        real ,      intent(INOUT) ::   pt(isd:ied  ,jsd:jed  ,npz)
        real ,      intent(INOUT) :: delp(isd:ied  ,jsd:jed  ,npz)
        real ,      intent(INOUT) ::    q(isd:ied  ,jsd:jed  ,npz, ncnst)
        
        real ,      intent(INOUT) :: phis(isd:ied  ,jsd:jed  )

        real ,      intent(INOUT) ::   ps(isd:ied  ,jsd:jed  )
        real ,      intent(INOUT) ::   pe(is-1:ie+1,npz+1,js-1:je+1)
        real ,      intent(INOUT) ::   pk(is:ie    ,js:je    ,npz+1)
        real ,      intent(INOUT) :: peln(is :ie   ,npz+1    ,js:je)
        real ,      intent(INOUT) ::  pkz(is:ie    ,js:je    ,npz  )
        
        real ,      intent(INOUT) ::   uc(isd:ied+1,jsd:jed  ,npz)
        real ,      intent(INOUT) ::   vc(isd:ied  ,jsd:jed+1,npz)
        real ,      intent(INOUT) ::   ua(isd:ied  ,jsd:jed  ,npz)
        real ,      intent(INOUT) ::   va(isd:ied  ,jsd:jed  ,npz)
        real ,      intent(inout) :: delz(is:ie,js:je,npz)
        real ,      intent(inout)   ::  ze0(is:ie,js:je,npz+1)
        
        real ,      intent(inout)    ::   ak(npz+1)
        real ,      intent(inout)    ::   bk(npz+1)
        
        integer,      intent(IN) :: npx, npy, npz
        integer,      intent(IN) :: ng, ncnst, nwat
        integer,      intent(IN) :: k_top
        integer,      intent(IN) :: ndims
        integer,      intent(IN) :: nregions
        
        real,         intent(IN) :: dry_mass
        logical,      intent(IN) :: mountain
        logical,      intent(IN) :: moist_phys
        logical,      intent(IN) :: hydrostatic, hybrid_z

        real, dimension(is:ie):: pm, qs
        real :: dist, r0, f0_const, prf, rgrav
        real :: ptmp, ze, zc, zm
        real :: t00, p00, xmax, xc, xx, yy, zz, pk0, pturb, ztop
        real :: ze1(npz+1)
         real:: dz1(npz)
        integer :: i, j, k, m, icenter, jcenter

        f0_const = 2.*omega*sin(deglat/180.*pi)
        f0(:,:) = f0_const
        fC(:,:) = f0_const

        q = 0.

        select case (test_case)
        case ( 1 )

           phis(:,:)=0.

           u (:,:,:)=10.
           v (:,:,:)=10.
           ua(:,:,:)=10.
           va(:,:,:)=10.
           uc(:,:,:)=10.
           vc(:,:,:)=10.
           pt(:,:,:)=1.
           delp(:,:,:)=0.
           
           do j=js,je
              if (j>0 .and. j<5) then
                 do i=is,ie
                    if (i>0 .and. i<5) then
                       delp(i,j,:)=1.
                    endif
                 enddo
              endif
           enddo
           call mpp_update_domains( delp, domain )

        case ( 2 )

           phis(:,:) = 0.

!          r0 = 5000.
           r0 = 5.*sqrt(dx_const**2 + dy_const**2)
           icenter = npx/2
           jcenter = npy/2
           do j=jsd,jed
              do i=isd,ied
                 dist=(i-icenter)*dx_const*(i-icenter)*dx_const   &
                       +(j-jcenter)*dy_const*(j-jcenter)*dy_const
                 dist=min(r0,sqrt(dist))
                 phis(i,j)=1500.*(1. - (dist/r0))
              enddo
           enddo

           u (:,:,:)=0.
           v (:,:,:)=0.
           ua(:,:,:)=0.
           va(:,:,:)=0.
           uc(:,:,:)=0.
           vc(:,:,:)=0.
           pt(:,:,:)=1.
           delp(:,:,:)=1500.

        case ( 14 )
!---------------------------
! Doubly periodic Aqua-plane
!---------------------------
           u(:,:,:) = 0.
           v(:,:,:) = 0.
           do j=jsd,jed
              do i=isd,ied
                 phis(i,j) = 0.
                   ps(i,j) = 1000.E2
              enddo
           enddo

           do k=1,npz
              do j=jsd,jed
                 do i=isd,ied
                    pt(i,j,k) = 250.  ! really cold start
                    delp(i,j,k) = ak(k+1)-ak(k) + ps(i,j)*(bk(k+1)-bk(k))
                 enddo
              enddo
           enddo

! *** Add Initial perturbation ***
           r0 = 20.*sqrt(dx_const**2 + dy_const**2)
!          icenter = npx/2
!          jcenter = npy/2
! Off center for spin up hurricanes
           icenter = npx/2 + 1
           jcenter = npy/2 + 1

           do j=js,je
              do i=is,ie
                 dist = (i-icenter)*dx_const*(i-icenter)*dx_const   &
                         +(j-jcenter)*dy_const*(j-jcenter)*dy_const
                 dist = min(r0,sqrt(dist))
                 do k=1,npz
                    prf = ak(k) + ps(i,j)*bk(k)
                    if ( prf > 100.E2 ) then
                         pt(i,j,k) = pt(i,j,k) + 50.*(1. - (dist/r0)) * prf/ps(i,j) 
                    endif
                 enddo
              enddo
           enddo

          call p_var(npz, is, ie, js, je, ptop, ptop_min, delp, delz, pt, ps,   &
                     pe, peln, pk, pkz, kappa, q, ng, ncnst, dry_mass, .false., .false., &
                     moist_phys, .true., k_top, nwat)

          q = 0.
         do k=3,npz
            do j=js,je
               do i=is,ie
                  pm(i) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
               enddo
               call qsmith(ie-is+1, 1, 1, pt(is:ie,j,k), pm, q(is:ie,j,k,1), qs)
               do i=is,ie
                  if ( pm(i)>100.E2 ) then
                       q(i,j,k,1) = 0.99*qs(i)
                  else
                       q(i,j,k,1) = 3.E-6
                  endif
               enddo
            enddo
         enddo

        case ( 15 )
!---------------------------
! Doubly periodic bubble
!---------------------------
           t00 = 250.

           u(:,:,:) = 0.
           v(:,:,:) = 0.
          pt(:,:,:) = t00
          q(:,:,:,:) = 1.E-6

          if ( .not. hydrostatic ) w(:,:,:) = 0.

           do j=jsd,jed
              do i=isd,ied
                 phis(i,j) = 0.
                   ps(i,j) = 1000.E2
              enddo
           enddo

           do k=1,npz
              do j=jsd,jed
                 do i=isd,ied
                    delp(i,j,k) = ak(k+1)-ak(k) + ps(i,j)*(bk(k+1)-bk(k))
                 enddo
              enddo
           enddo


           do k=1,npz
              do j=jsd,jed
                 do i=isd,ied
                         ptmp = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
!                   pt(i,j,k) = t00
                 enddo
              enddo
           enddo

          call p_var(npz, is, ie, js, je, ptop, ptop_min, delp, delz, pt, ps,   &
                     pe, peln, pk, pkz, kappa, q, ng, ncnst, dry_mass, .false., .false., &
                     moist_phys, .false., k_top, nwat)

! *** Add Initial perturbation ***
           r0 = 5.*max(dx_const, dy_const)
           zc = 0.5e3         ! center of bubble  from surface
           icenter = npx/2
           jcenter = npy/2

           do j=js,je
              do i=is,ie
                 ze = 0.
                 do k=npz,1,-1
                    zm = ze - 0.5*delz(i,j,k)   ! layer center
                    ze = ze - delz(i,j,k)
                    dist = ((i-icenter)*dx_const)**2 + ((j-jcenter)*dy_const)**2 +  &
                           (zm-zc)**2
                    dist = sqrt(dist)
                    if ( dist <= r0 ) then
                         pt(i,j,k) = pt(i,j,k) + 5.*(1.-dist/r0)
                    endif
                 enddo
              enddo
           enddo

        case ( 16 )
!------------------------------------
! Non-hydrostatic 3D density current:
!------------------------------------
           phis = 0.
           u = 0.
           v = 0.
           w = 0.
          t00 = 300.
          p00 = 1.E5
          pk0 = p00**kappa
! Set up vertical coordinare with constant del-z spacing:
! Control: npz=64;  dx = 100 m; dt = 1; n_split=10
          ztop = 6.4E3
         ze1(    1) = ztop
         ze1(npz+1) = 0.
         do k=npz,2,-1
            ze1(k) = ze1(k+1) + ztop/real(npz)
         enddo

          do j=js,je
             do i=is,ie
                ps(i,j) = p00
                pe(i,npz+1,j) = p00
                pk(i,j,npz+1) = pk0
             enddo
          enddo

          do k=npz,1,-1
             do j=js,je
                do i=is,ie
                   delz(i,j,k) = ze1(k+1) - ze1(k)
                     pk(i,j,k) = pk(i,j,k+1) + grav*delz(i,j,k)/(cp_air*t00)*pk0
                     pe(i,k,j) = pk(i,j,k)**(1./kappa) 
                enddo
             enddo
          enddo

          ptop = pe(is,1,js)
          if ( gid==masterproc ) write(*,*) 'Density curent testcase: model top (mb)=', ptop/100.

          do k=1,npz+1
             do j=js,je
                do i=is,ie
                   peln(i,k,j) = log(pe(i,k,j)) 
                    ze0(i,j,k) = ze1(k)
                enddo
             enddo
          enddo

          do k=1,npz
             do j=js,je
                do i=is,ie
                   pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
                  delp(i,j,k) =  pe(i,k+1,j)-pe(i,k,j)  
                    pt(i,j,k) = t00/pk0   ! potential temp
                enddo
             enddo
          enddo

          pturb = 15.
           xmax = 51.2E3 
             xc = xmax / 2.

         do k=1,npz
            zz = (0.5*(ze1(k)+ze1(k+1))-3.E3) / 2.E3
            do j=js,je
               do i=is,ie
! Impose perturbation in potential temperature: pturb
                  xx = (dx_const * (0.5+real(i-1)) - xc) / 4.E3 
                  yy = (dy_const * (0.5+real(j-1)) - xc) / 4.E3
                  dist = sqrt( xx**2 + yy**2 + zz**2 )
                  if ( dist<=1. ) then
                       pt(i,j,k) = pt(i,j,k) - pturb/pkz(i,j,k)*(cos(pi*dist)+1.)/2. 
                  endif
! Transform back to temperature:
                  pt(i,j,k) = pt(i,j,k) * pkz(i,j,k)
               enddo
            enddo
          enddo

        case ( 101 )

! IC for LES
         t00 = 250.      ! constant temp
         p00 = 1.E5
         pk0 = p00**kappa

         phis = 0.
         u = 0.
         v = 0.
         w = 0.
         pt(:,:,:) = t00
         q(:,:,:,1) = 0.

         if (.not.hybrid_z) call mpp_error(FATAL, 'hybrid_z must be .TRUE.')

          rgrav = 1./ grav

          if ( npz/=101)  then
              call mpp_error(FATAL, 'npz must be == 101 ')
          else
              call compute_dz_L101( npz, ztop, dz1 )
          endif

          call set_hybrid_z(is, ie, js, je, ng, npz, ztop, dz1, rgrav,  &
                            phis, ze0, delz)

          do j=js,je
             do i=is,ie
                ps(i,j) = p00
                pe(i,npz+1,j) = p00
                pk(i,j,npz+1) = pk0
                peln(i,npz+1,j) = log(p00)
             enddo
          enddo

          do k=npz,1,-1
             do j=js,je
                do i=is,ie
                   peln(i,k,j) = peln(i,k+1,j) + grav*delz(i,j,k)/(rdgas*t00)
                     pe(i,k,j) = exp(peln(i,k,j))
                     pk(i,j,k) = pe(i,k,j)**kappa
                enddo
             enddo
          enddo


! Set up fake "sigma" coordinate 
          call make_eta_level(npz, pe, area, ks, ak, bk)

          if ( gid==masterproc ) write(*,*) 'LES testcase: computed model top (mb)=', ptop/100.

          do k=1,npz
             do j=js,je
                do i=is,ie
                   pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
                  delp(i,j,k) =  pe(i,k+1,j)-pe(i,k,j)  
                enddo
             enddo
          enddo

         do k=1,npz
            do j=js,je
               do i=is,ie
                  pm(i) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
               enddo
               call qsmith(ie-is+1, 1, 1, pt(is:ie,j,k), pm, q(is:ie,j,k,1), qs)
               do i=is,ie
                  if ( pm(i) > 100.E2 ) then
                       q(i,j,k,1) = 0.9*qs(i)
                  else
                       q(i,j,k,1) = 2.E-6
                  endif
               enddo
            enddo
         enddo

! *** Add perturbation ***
           r0 = 1.0e3         ! radius (m)
           zc = 1.0e3         ! center of bubble 
           icenter = npx/2
           jcenter = npy/2

           do k=1,npz
              do j=js,je
                 do i=is,ie
                    zm = 0.5*(ze0(i,j,k)+ze0(i,j,k+1))
                    dist = ((i-icenter)*dx_const)**2 + ((j-jcenter)*dy_const)**2 + (zm-zc)**2
                    dist = sqrt(dist)
                    if ( dist <= r0 ) then
                         pt(i,j,k) = pt(i,j,k) + 2.0*(1.-dist/r0)
                    endif
                 enddo
              enddo
           enddo

        end select

      end subroutine init_double_periodic

      subroutine init_latlon(u,v,pt,delp,q,phis, ps,pe,peln,pk,pkz,  uc,vc, ua,va, ak, bk,  &
                             npx, npy, npz, ng, ncnst, k_top, ndims, nregions, dry_mass,    &
                             mountain, moist_phys, hybrid_z, delz, ze0)

        real ,      intent(INOUT) ::    u(isd:ied  ,jsd:jed+1,npz)
        real ,      intent(INOUT) ::    v(isd:ied+1,jsd:jed  ,npz)
        real ,      intent(INOUT) ::   pt(isd:ied  ,jsd:jed  ,npz)
        real ,      intent(INOUT) :: delp(isd:ied  ,jsd:jed  ,npz)
        real ,      intent(INOUT) ::    q(isd:ied  ,jsd:jed  ,npz, ncnst)
        
        real ,      intent(INOUT) :: phis(isd:ied  ,jsd:jed  )

        real ,      intent(INOUT) ::   ps(isd:ied  ,jsd:jed  )
        real ,      intent(INOUT) ::   pe(is-1:ie+1,npz+1,js-1:je+1)
        real ,      intent(INOUT) ::   pk(is:ie    ,js:je    ,npz+1)
        real ,      intent(INOUT) :: peln(is :ie   ,npz+1    ,js:je)
        real ,      intent(INOUT) ::  pkz(is:ie    ,js:je    ,npz  )
        
        real ,      intent(INOUT) ::   uc(isd:ied+1,jsd:jed  ,npz)
        real ,      intent(INOUT) ::   vc(isd:ied  ,jsd:jed+1,npz)
        real ,      intent(INOUT) ::   ua(isd:ied  ,jsd:jed  ,npz)
        real ,      intent(INOUT) ::   va(isd:ied  ,jsd:jed  ,npz)
        real ,      intent(inout) :: delz(is:ie,js:je,npz)
        real ,      intent(inout)   ::  ze0(is:ie,js:je,npz+1)
        
        real ,      intent(IN)    ::   ak(npz+1)
        real ,      intent(IN)    ::   bk(npz+1)
        
        integer,      intent(IN) :: npx, npy, npz
        integer,      intent(IN) :: ng, ncnst
        integer,      intent(IN) :: k_top
        integer,      intent(IN) :: ndims
        integer,      intent(IN) :: nregions
        
        real,         intent(IN) :: dry_mass
        logical,      intent(IN) :: mountain
        logical,      intent(IN) :: moist_phys
        logical,      intent(IN) :: hybrid_z

        real    :: p1(2), p2(2), r, r0
        integer :: i,j

        do j=jsd,jed+1
           do i=isd,ied+1
              fc(i,j) = 2.*omega*( -cos(grid(i,j,1))*cos(grid(i,j,2))*sin(alpha)  &
                                   +sin(grid(i,j,2))*cos(alpha) )
           enddo
        enddo
        do j=jsd,jed
           do i=isd,ied
              f0(i,j) = 2.*omega*( -cos(agrid(i,j,1))*cos(agrid(i,j,2))*sin(alpha)  &
                                   +sin(agrid(i,j,2))*cos(alpha) )
           enddo
        enddo

        select case (test_case)
        case ( 1 )

         Ubar = (2.0*pi*radius)/(12.0*86400.0)
         phis = 0.0
         r0 = radius/3. !RADIUS radius/3.
!!$         p1(1) = 0.
         p1(1) = pi/2. + pi_shift
         p1(2) = 0.
         do j=jsd,jed
            do i=isd,ied
               p2(1) = agrid(i,j,1)
               p2(2) = agrid(i,j,2)
               r = great_circle_dist( p1, p2, radius )
               if (r < r0) then
                  delp(i,j,1) = phis(i,j) + 0.5*(1.0+cos(PI*r/r0))
               else
                  delp(i,j,1) = phis(i,j)
               endif
            enddo
         enddo
         call init_latlon_winds(UBar, u, v, ua, va, uc, vc, 1)


!!$           phis(:,:)=0.
!!$
!!$           u (:,:,:)=10.
!!$           v (:,:,:)=10.
!!$           ua(:,:,:)=10.
!!$           va(:,:,:)=10.
!!$           uc(:,:,:)=10.
!!$           vc(:,:,:)=10.
!!$           pt(:,:,:)=1.
!!$           delp(:,:,:)=0.
!!$           
!!$           do j=js,je
!!$              if (j>10 .and. j<15) then
!!$                 do i=is,ie
!!$                    if (i>10 .and. i<15) then
!!$                       delp(i,j,:)=1.
!!$                    endif
!!$                 enddo
!!$              endif
!!$           enddo
!!$           call mpp_update_domains( delp, domain )

        end select

      end subroutine init_latlon

      subroutine init_latlon_winds(UBar, u, v, ua, va, uc, vc, defOnGrid)

        ! defOnGrid = -1:null_op, 0:All-Grids, 1:C-Grid, 2:D-Grid, 3:A-Grid, 4:A-Grid then Rotate, 5:D-Grid with unit vectors then Rotate

        real,    intent(INOUT) :: UBar
        real,    intent(INOUT) ::  u(isd:ied  ,jsd:jed+1)
        real,    intent(INOUT) ::  v(isd:ied+1,jsd:jed  )
        real,    intent(INOUT) :: uc(isd:ied+1,jsd:jed  )
        real,    intent(INOUT) :: vc(isd:ied  ,jsd:jed+1)
        real,    intent(INOUT) :: ua(isd:ied  ,jsd:jed  )
        real,    intent(INOUT) :: va(isd:ied  ,jsd:jed  )
        integer, intent(IN)    :: defOnGrid

        real   :: p1(2),p2(2),p3(2),p4(2), pt(2)
        real :: e1(3), e2(3), ex(3), ey(3)

        real   :: dist, r, r0 
        integer :: i,j,k,n
        real :: utmp, vtmp

        real :: psi_b(isd:ied+1,jsd:jed+1), psi(isd:ied,jsd:jed), psi1, psi2 

        psi(:,:) = 1.e25
        psi_b(:,:) = 1.e25
        do j=jsd,jed
           do i=isd,ied
              psi(i,j) = (-1.0 * Ubar * radius *( sin(agrid(i,j,2))                  *cos(alpha) - &
                                                  cos(agrid(i,j,1))*cos(agrid(i,j,2))*sin(alpha) ) )
           enddo
        enddo
        do j=jsd,jed+1
           do i=isd,ied+1
              psi_b(i,j) = (-1.0 * Ubar * radius *( sin(grid(i,j,2))                 *cos(alpha) - &
                                                    cos(grid(i,j,1))*cos(grid(i,j,2))*sin(alpha) ) )
           enddo
        enddo
        
        if ( defOnGrid == 1 ) then
           do j=jsd,jed+1
              do i=isd,ied
                 dist = dx(i,j)
                 vc(i,j) = (psi_b(i+1,j)-psi_b(i,j))/dist
                 if (dist==0) vc(i,j) = 0.
              enddo
           enddo
           do j=jsd,jed
              do i=isd,ied+1
                 dist = dy(i,j)
                 uc(i,j) = -1.0*(psi_b(i,j+1)-psi_b(i,j))/dist
                 if (dist==0) uc(i,j) = 0.
              enddo
           enddo

           
           do j=js,je
              do i=is,ie+1
                 dist = dxc(i,j)
                 v(i,j) = (psi(i,j)-psi(i-1,j))/dist
                 if (dist==0) v(i,j) = 0.            
              enddo
           enddo
           do j=js,je+1
              do i=is,ie
                 dist = dyc(i,j)
                 u(i,j) = -1.0*(psi(i,j)-psi(i,j-1))/dist
                 if (dist==0) u(i,j) = 0. 
              enddo
           enddo
        endif
     
      end subroutine init_latlon_winds
      
end module test_cases_mod
