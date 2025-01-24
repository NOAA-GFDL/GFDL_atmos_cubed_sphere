!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

     module fv_ideal_mod

#ifdef OVERLOAD_R4
      use constantsR4_mod,   only: pi=>pi_8, omega, grav, kappa, rdgas, cp_air, rvgas
#else
      use constants_mod,     only: pi=>pi_8, omega, grav, kappa, rdgas, cp_air, rvgas
#endif
      use init_hydro_mod,    only: p_var
      use fv_mp_mod,         only: is_master
      use fv_grid_utils_mod, only: ptop_min
      use mpp_mod,           only: mpp_error, FATAL, stdlog, input_nml_file
      use fms_mod,           only: check_nml_error     
      use mpp_domains_mod,   only: domain2d
      use fv_arrays_mod,         only: fv_grid_type, fv_flags_type, fv_grid_bounds_type, R_GRID
      implicit none
      private

      integer :: sphum, theta_d
      integer, parameter :: max_bub=100
      real(kind=R_GRID), parameter :: one = 1.d0
      integer :: test_case = 11
      logical :: bubble_do = .false.
      logical :: do_rand_perts = .false.
      real    :: alpha = 0.0
      integer :: Nsolitons = 1, n_bub = 1
      real    :: soliton_size = 750.e3, soliton_Umax = 50.
      integer :: t_profile = 0, q_profile = 0, ws_profile = 0
      integer :: do_coriolis = 0, bubble_type = 0
      real    :: bubble_t = 2., bubble_q = 0., bubble_rad_x = 10.0E3
      real    :: umove = 0.0, vmove = 0.0
      real    :: p00_in = 1.e5
      real    :: bubble_rad_y = 10.0E3, bubble_zc = 1.4E3
      real    :: iso_t = 300., adi_th = 300., us0 = 30.
      real,dimension(max_bub)    :: icenters, jcenters

     public :: fv_init_ideal
     public :: read_namelist_fv_ideal
     contains

     subroutine fv_init_ideal(u,v,w,pt,delp,q,phis, ps,pe,peln,pk,pkz,  uc,vc, ua,va, ak, bk,  &
                              gridstruct, flagstruct, npx, npy, npz, ng, ncnst, nwat, ndims, nregions, dry_mass, &
                              mountain, moist_phys, hydrostatic, hybrid_z, delz, ze0, ks, ptop, domain_in, tile_in, bd)

        type(fv_grid_bounds_type), intent(IN) :: bd
        real ,      intent(INOUT) ::    u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)
        real ,      intent(INOUT) ::    v(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)
        real ,      intent(INOUT) ::    w(bd%isd:  ,bd%jsd:  ,1:)
        real ,      intent(INOUT) ::   pt(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
        real ,      intent(INOUT) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
        real ,      intent(INOUT) ::    q(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst)

        real ,      intent(INOUT) :: phis(bd%isd:bd%ied  ,bd%jsd:bd%jed  )

        real ,      intent(INOUT) ::   ps(bd%isd:bd%ied  ,bd%jsd:bd%jed  )
        real ,      intent(INOUT) ::   pe(bd%is-1:bd%ie+1,npz+1,bd%js-1:bd%je+1)
        real ,      intent(INOUT) ::   pk(bd%is:bd%ie    ,bd%js:bd%je    ,npz+1)
        real ,      intent(INOUT) :: peln(bd%is :bd%ie   ,npz+1    ,bd%js:bd%je)
        real ,      intent(INOUT) ::  pkz(bd%is:bd%ie    ,bd%js:bd%je    ,npz  )
        real ,      intent(INOUT) ::   uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)
        real ,      intent(INOUT) ::   vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)
        real ,      intent(INOUT) ::   ua(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
        real ,      intent(INOUT) ::   va(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
        real ,      intent(inout) :: delz(bd%is:,bd%js:,1:)
        real ,      intent(inout)   ::  ze0(bd%is:,bd%js:,1:)

        real ,      intent(inout)    ::   ak(npz+1)
        real ,      intent(inout)    ::   bk(npz+1)

        integer,      intent(IN) :: npx, npy, npz
        integer,      intent(IN) :: ng, ncnst, nwat
        integer,      intent(IN) :: ndims
        integer,      intent(IN) :: nregions

        real,         intent(IN) :: dry_mass
        logical,      intent(IN) :: mountain
        logical,      intent(IN) :: moist_phys
        logical,      intent(IN) :: hydrostatic, hybrid_z
        integer,      intent(INOUT) :: ks
        integer,      intent(INOUT), target :: tile_in
        real,         intent(INOUT) :: ptop

        type(domain2d), intent(IN), target :: domain_in

        type(fv_grid_type), intent(IN), target :: gridstruct
        type(fv_flags_type), intent(IN), target :: flagstruct

        integer, parameter :: nl_max = 2500
        real, dimension(nl_max) ::  z_snd, p_snd, t_snd, rho_snd, u_snd, v_snd, qv_snd
        real, dimension(bd%is:bd%ie):: pm, qs
        real, dimension(1:npz):: pk1, pe1, ts1, qs1, dummy
        !real :: us0 = 30.
        real :: dist,r0, f0_const, prf, rgrav, xrad, yrad, zrad,RAD
        real :: xradbub, yradbub
        real :: ptmp, ze, zc, zm, utmp, vtmp
        real :: t00, p00, xmax, xc, xx, yy, pk0, pturb, ztop
        real :: ze1(npz+1)
        real:: dz1(npz)
        real:: zvir, rand1, rand2
        real:: amplitude = 0.2
        integer :: i, j, k, m, icenter, jcenter, nl_snd

        real, pointer, dimension(:,:,:)   :: agrid, grid
        real(kind=R_GRID), pointer, dimension(:,:)     :: area
        real, pointer, dimension(:,:)     :: rarea, fC, f0
        real, pointer, dimension(:,:,:)   :: ee1, ee2, en1, en2
        real, pointer, dimension(:,:,:,:) :: ew, es
        real, pointer, dimension(:,:)     :: dx,dy, dxa,dya, rdxa, rdya, dxc,dyc

        logical, pointer :: cubed_sphere, latlon

        type(domain2d), pointer :: domain
        integer, pointer :: tile

        logical, pointer :: have_south_pole, have_north_pole

        integer, pointer :: ntiles_g
        real,    pointer :: acapN, acapS, globalarea

        real(kind=R_GRID), pointer :: dx_const, dy_const

        integer :: is,  ie,  js,  je
        integer :: isd, ied, jsd, jed

        integer :: b

        is  = bd%is
        ie  = bd%ie
        js  = bd%js
        je  = bd%je
        isd = bd%isd
        ied = bd%ied
        jsd = bd%jsd
        jed = bd%jed

        agrid => gridstruct%agrid
        grid  => gridstruct%grid

        area => gridstruct%area_64

        dx      => gridstruct%dx
        dy      => gridstruct%dy
        dxa     => gridstruct%dxa
        dya     => gridstruct%dya
        rdxa    => gridstruct%rdxa
        rdya    => gridstruct%rdya
        dxc     => gridstruct%dxc
        dyc     => gridstruct%dyc

        fC    => gridstruct%fC
        f0    => gridstruct%f0

        !These are frequently used and so have pointers set up for them
        dx_const => flagstruct%dx_const
        dy_const => flagstruct%dy_const

        domain => domain_in
        tile => tile_in

        have_south_pole               => gridstruct%have_south_pole
        have_north_pole               => gridstruct%have_north_pole

        ntiles_g                      => gridstruct%ntiles_g
        acapN                         => gridstruct%acapN
        acapS                         => gridstruct%acapS
        globalarea                    => gridstruct%globalarea

        if (do_coriolis>=1) then
           f0_const = 2.*omega*sin(flagstruct%deglat/180.*pi)
        else
           f0_const = 0.0
        endif

        f0(:,:) = f0_const
        fC(:,:) = f0_const

        q = 0.

        zvir = rvgas/rdgas - 1.
        p00 = p00_in ! 1000.E2
        print*, do_coriolis, t_profile, q_profile, ws_profile
        if (t_profile == -1 .or. q_profile == -1 .or. ws_profile == -1) then
          call get_sounding( z_snd, p_snd, t_snd, rho_snd, u_snd, v_snd, qv_snd, nl_max, nl_snd, p00)
          IF ( is_master() ) write(*,*) 'Using p00 from sounding file: p00 = ',p00
        endif

   ! Set up pressure arrays
        ps(:,:) = p00
        phis(:,:) = 0.
        do j=js,je
           do i=is,ie
                pk(i,j,1) = ptop**kappa
                pe(i,1,j) = ptop
              peln(i,1,j) = log(ptop)
           enddo
        enddo

        do k=1,npz
           do j=js,je
              do i=is,ie
                 delp(i,j,k) = ak(k+1)-ak(k) + ps(i,j)*(bk(k+1)-bk(k))
                 pe(i,k+1,j) = ak(k+1) + ps(i,j)*bk(k+1)
                 peln(i,k+1,j) = log(pe(i,k+1,j))
                 pk(i,j,k+1) = exp( kappa*peln(i,k+1,j) )
              enddo
           enddo
        enddo

        i = is
        j = js
        do k=1,npz
           pk1(k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
           pe1(k) = (pe(i,k+1,j)-pe(i,k,j))/(peln(i,k+1,j)-peln(i,k,j))
        enddo

   !Set model thermodynamic profile based on namelist inputs
        ! Read profile from sounding
        if (t_profile == -1) then
          do k = 1,npz
           ts1(k) =  interp_log( t_snd, p_snd, pe1(k), nl_max, nl_snd  )
          enddo
        elseif ( t_profile == 0 ) then
        !Generate GFDL's supercell sounding
          call SuperCell_Sounding(npz, p00, pk1, ts1, qs1)
        elseif (t_profile == 1 ) then
          !adiabatic
          print*, "kappa, adi_th = ", kappa, adi_th
          do k=1,npz
            ts1(k) = adi_th * ( pe1(k) / 1E5) ** kappa
          enddo
        elseif (t_profile == 2) then
          !isothermal
          ts1(:) = iso_t
        else
          call mpp_error(FATAL, " t_profile ", t_profile ,"  not defined" )
        endif

        if (q_profile == -1) then
          ! Read profile from sounding
          do k = 1,npz
           qs1(k) =  interp_log( qv_snd, p_snd, pe1(k), nl_max,nl_snd  )
          enddo
        elseif ( q_profile==0 ) then
          if (t_profile == 0) then
            ! qs1 already computed prior, move along
          else
            ! Generate GFDL's supercell sounding
            call SuperCell_Sounding(npz, p00, pk1, dummy, qs1)
          endif
        elseif (q_profile==1 ) then
          ! dry environment
          qs1(:) = 1E-9
        else
          call mpp_error(FATAL, " q_profile ", q_profile ,"  not defined" )
        endif

        ! Compute delz from ts1 and qs1
        w(:,:,:) = 0.
        q(:,:,:,:) = 0.

        do k=1,npz
           do j=js,je
              do i=is,ie
                 pt(i,j,k)   = ts1(k)
                  q(i,j,k,1) = qs1(k)
                 delz(i,j,k)=rdgas/grav*ts1(k)*(1.+zvir*qs1(k))*(peln(i,k,j)-peln(i,k+1,j))
                enddo
             enddo
          enddo

        ze1(npz+1) = 0.
        do k=npz,1,-1
           ze1(k) = ze1(k+1) - delz(is,js,k)
        enddo

    ! Set up model winds
        !Read winds from sounding
        if (ws_profile == -1) then
          do k = 1,npz
           zm = 0.5*(ze1(k)+ze1(k+1))
           utmp =  interp_lin( u_snd, z_snd, zm, nl_max, nl_snd  )
           vtmp = interp_lin( v_snd, z_snd,zm, nl_max, nl_snd  )
            do j=js,je+1
              do i=is,ie
                 u(i,j,k) = utmp
             enddo
           enddo
            do j=js,je
              do i=is,ie+1
                v(i,j,k) = vtmp
              enddo
            enddo
          enddo

        elseif ( ws_profile==0 ) then
        ! Quarter-circle hodograph (Harris approximation)
          do k=1,npz
           zm = 0.5*(ze1(k)+ze1(k+1))
           if ( zm .le. 2.e3 ) then
               utmp = 8.*(1.-cos(pi*zm/4.e3))
               vtmp = 8.*sin(pi*zm/4.e3)
           elseif (zm .le. 6.e3 ) then
               utmp = 8. + (us0-8.)*(zm-2.e3)/4.e3
               vtmp = 8.
           else
               utmp = us0
               vtmp = 8.
           endif
! u-wind
           do j=js,je+1
              do i=is,ie
                 u(i,j,k) = utmp - 8.
             enddo
           enddo
! v-wind
           do j=js,je
              do i=is,ie+1
                 v(i,j,k) = vtmp - 4.
             enddo
           enddo
          enddo
        elseif (ws_profile==1 ) then
        ! Unidirectional WK shear
          v(:,:,:) = 0.
          do k=1,npz
            zm = 0.5*(ze1(k)+ze1(k+1))
            if ( zm .le. 6.e3 ) then
              u(:,:,k) = us0 * tanh(zm/3.e3)
            else
              u(:,:,k) = us0
            endif
          enddo
        elseif (ws_profile==2 ) then
        ! constant u, v=0
          u(:,:,:) = us0
          v(:,:,:) = 0.
        elseif (ws_profile==3 ) then
        ! constant v, u=0
          u(:,:,:) = 0.
          v(:,:,:) = us0
        elseif (ws_profile==4 ) then
          u(:,:,:) = 0.
          v(:,:,:) = 0.
        elseif (ws_profile==5) then
          ! Linear WK shear below 2.5km
          v(:,:,:) = 0.
          do k=1,npz
            zm = 0.5*(ze1(k)+ze1(k+1))
            if ( zm .le. 2.5e3 ) then
              u(:,:,k) = us0 * tanh(zm/3.e3)
            else
              u(:,:,k) = us0
            endif
          enddo
        else
         call mpp_error(FATAL, " ws_profile ", ws_profile ,"  not defined" )
        endif

        IF ( is_master() ) THEN
          write(*,*) 'Final sounding: k, z, t, q, u, v, p'
          do k = 1,npz
            write(*,*) k,ze1(k),ts1(k),qs1(k),u(is,js,k),v(is,js,k),pe1(k)
          enddo
        ENDIF

        call p_var(npz, is, ie, js, je, ptop, ptop_min, delp, delz, pt, ps,&
                   pe, peln, pk, pkz, kappa, q, ng, ncnst, area, dry_mass,.false.,.false., &
                  .true., hydrostatic, nwat, domain, flagstruct%adiabatic)

        ! Add in (or don't) bubble(s) to initialize storms
        if (bubble_type > 0) then
        if (is_master()) print*, "ADDING BUBBLE"
! *** Add Initial perturbation ***
          pturb = bubble_t
          xradbub = bubble_rad_x
          yradbub = bubble_rad_y
          zc = bubble_zc     ! center of bubble from surface
          if (bubble_type == 1) then
            ! One bubble in domain center
            icenter = (npx-1)/2 + 1
            jcenter = (npy-1)/2 + 1
            n_bub = 1
          elseif (bubble_type == 2) then
            !Line of centered N-S bubbles for squall line
            !n_bub = floor(float(npy)*dy_const/r0)
            if ( is_master() )  print*, "initializing ", n_bub , " bubbles"
            icenter = 0
          elseif (bubble_type == 3) then
            ! User entry of i/j bubble locations
            if ( is_master() ) print*, "initializing ", n_bub , " bubbles"
            if ( is_master() ) print*, "at locations i = ", icenters(1:n_bub), "j = ", jcenters(1:n_bub)
          endif
          do j = js, je
          do i = is, ie
          do k=1, npz
            do b = 1, n_bub
              call random_number(rand1)
              call random_number(rand2)
              if (bubble_type == 2) then
                jcenter=((npy-1)/2+1)-((n_bub+1)/2-b)*30000.0/dy_const
              elseif (bubble_type == 3) then
                icenter = icenters(b)
                jcenter = jcenters(b)
              endif
              zm = 0.5*(ze1(k)+ze1(k+1))
              yrad = dy_const*float(j-jcenter)/yradbub
              xrad = dx_const*float(i-icenter)/xradbub

              RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
              IF(RAD <= 1.) THEN
                 if (do_rand_perts) then
                    ! Add in random, small amplitude perturbations to bubble thermodynamic state
                    pt(i,j,k) = pt(i,j,k) + pturb*COS(.5*pi*RAD)**2 + 0.2 * (2.0*rand1-1.0)
                    q(i,j,k,1) = q(i,j,k,1) + bubble_q *COS(.5*pi*RAD)**2 + 1.0E-7 *(2.0*rand2-1.0)
                 else
                    pt(i,j,k) = pt(i,j,k) + pturb*COS(.5*pi*RAD)**2
                    q(i,j,k,1) = q(i,j,k,1) + bubble_q *COS(.5*pi*RAD)**2
                 endif
              ENDIF
             enddo !nbub
           enddo!npz
           enddo!i
           enddo!j
        endif !bubbletype

        uc(isd:ied,:,:) =  u(:,jsd:jed,:)
        uc(ied+1,:,:) = u(ied,jsd:jed,:)
        ua(:,:,:) = u(:,jsd:jed,:)

        vc(:,jsd:jed,:) = v(isd:ied,:,:)
        vc(:,jed+1,:) = v(isd:ied,jed,:)
        va(:,:,:) = v(isd:ied,:,:)

        nullify(grid)
        nullify(agrid)
   
        nullify(area)
   
        nullify(fC)
        nullify(f0)
    
        nullify(ee1)
        nullify(ee2)
        nullify(ew)
        nullify(es)
        nullify(en1)
        nullify(en2)
   
        nullify(dx)
        nullify(dy)
        nullify(dxa)
        nullify(dya)
        nullify(rdxa)
        nullify(rdya)
        nullify(dxc)
        nullify(dyc)
   
        nullify(dx_const)
        nullify(dy_const)
   
        nullify(domain)
        nullify(tile)
   
        nullify(have_south_pole)
        nullify(have_north_pole)
   
        nullify(ntiles_g)
        nullify(acapN)
        nullify(globalarea)

     end subroutine fv_init_ideal

     subroutine read_namelist_fv_ideal(nml_filename)

        character(*), intent(IN) :: nml_filename
        integer :: ierr, f_unit, unit, ios
        namelist /test_case_nml/bubble_do, &
                                t_profile, q_profile, ws_profile, bubble_t, bubble_q,  &
                                bubble_zc, do_coriolis, iso_t, adi_th, us0, bubble_type,n_bub, &
                                icenters,jcenters, bubble_rad_x, bubble_rad_y, do_rand_perts, &
                                p00_in, umove, vmove

#include<file_version.h>

        unit = stdlog()

        ! Make alpha = 0 the default:
        alpha = 0.
        bubble_do = .false.
        test_case = 11   ! (USGS terrain)

        ! Read Test_Case namelist
        read (input_nml_file,test_case_nml,iostat=ios)
        ierr = check_nml_error(ios,'test_case_nml')
        write(unit, nml=test_case_nml)


      end subroutine read_namelist_fv_ideal 

     subroutine get_sounding( zk, p, t, rho, u, v, qv, nl_max, nl_in, p00 )
      implicit none

      integer nl_max, nl_in
      real p00
      real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
           u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max), t(nl_max)

      integer n
      parameter(n=1000)
      logical debug
      parameter( debug = .true.)
      character*256 message

! input sounding data

      real p_surf, th_surf, qv_surf
      real pi_surf, pi(n)
      real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

! diagnostics

      real rho_surf, p_input(n), rho_input(n)
      real pm_input(n)  !  this are for full moist sounding

! local data

      real r, g,cp
      parameter (r = rdgas)
      parameter (g = grav)
      parameter (cp = cp_air)
      real, parameter :: p1000mb = 100000.0, cvpm = -718./cp_air
      real, parameter :: rvovrd = rvgas/rdgas
      integer k, it, nl, nl_file, istat
      real qvf, qvf1, dz
      character*256 line

!  first, read the sounding
     print*,"OPEN SOUNDING FILE: ", trim('input_sounding')
     open(14, file=trim('input_sounding'), form='formatted', iostat=istat)
     if (istat /= 0) then
       call mpp_error(FATAL,"ERROR OPENING VARIABLE MAPPING FILE")
     endif
  
     nl = 0
     nl_file = 0
  
     !Loop over lines of file to count the number of levels
     do
       read(14, '(A)', iostat=istat) line
       if (istat/=0) exit
       if ( trim(line) .eq. '' ) cycle
       nl_file = nl_file + 1
     enddo
     if ( nl_file == 0) call mpp_error(FATAL,"VARMAP FILE IS EMPTY.")
     nl = nl_file -1
  
     nl_in = nl
     if(nl_in .gt. nl_max ) then
       print*, ' too many levels for input arrays ',nl_in,nl_max
       call mpp_error (FATAL, 'Too many levels for input arrays ' )
      end if
  
     rewind(14)
     read(14,*,iostat=istat) p_surf, th_surf, qv_surf
      do k = 2,nl_file
        read(14, *, iostat=istat) h_input(nl_in-k+2), th_input(nl_in-k+2), qv_input(nl_in-k+2), u_input(nl_in-k+2), v_input(nl_in-k+2)
       if (istat /= 0) call mpp_error(FATAL,"READING VARIABLE MAPPING FILE")
     enddo
     close(14)
  
     !  compute diagnostics,
  !  first, convert qv(g/kg) to qv(g/g)
  
      if (qv_surf > 0.1) THEN
        do k=1,nl
          qv_input(k) = 0.001*qv_input(k)
        enddo
      endif
  
      do k=1,nl
         IF ( qv_input(k) < 1E-9 ) THEN
           write(*,*) 'Warning Input sounding has qv = 0, resetting to 1E-9 kg/kg'
           qv_input(k) = 1E-9 ! Max(0.00001,qv_input(k))
         ENDIF
      enddo
  
      IF ( p_surf < 2000. ) p_surf = 100.*p_surf  ! convert to pascals
      p00 = p_surf
      qvf = 1. + rvovrd*qv_input(1)
      rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
      pi_surf = (p_surf/p1000mb)**(r/cp)
  
  
  
  !  integrate moist sounding hydrostatically, starting from the
  !  specified surface pressure
  !  -> first, integrate from surface to lowest level
  
            qvf = 1. + rvovrd*qv_input(nl)
            qvf1 = 1. + qv_input(nl)
            rho_input(nl) = rho_surf
            dz = h_input(nl)
            do it=1,10
              pm_input(nl) = p_surf &
                      - 0.5*dz*(rho_surf+rho_input(nl))*g*qvf1
              rho_input(nl) = 1./((r/p1000mb)*th_input(nl)*qvf*((pm_input(nl)/p1000mb)**cvpm))
            enddo
            do k=nl-1,1,-1
              rho_input(k) = rho_input(k+1)
              dz = h_input(k)-h_input(k+1)
              qvf1 = 0.5*(2.+(qv_input(k+1)+qv_input(k)))
              qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here
  
              do it=1,10
                pm_input(k) = pm_input(k+1) &
                        - 0.5*dz*(rho_input(k)+rho_input(k+1))*g*qvf1
                IF(pm_input(k) .LE. 0. )THEN
                  print*, "Integrated pressure has gone negative - toocold for chosen height"
                  WRITE(message,*)'k,pm_input(k),h_input(k),th_input(k) =',k,pm_input(k),h_input(k),th_input(k)
                  CALL mpp_error (FATAL, message )
                ENDIF
                rho_input(k) =1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
              enddo
            enddo
  
         do k=1,nl
  
            zk(k) = h_input(k)
            p(k) = pm_input(k)
            t(k) = th_input(k) * (pm_input(k) / p1000mb)**(rdgas/cp_air)
            u(k) = u_input(k) - umove
            v(k) = v_input(k) - vmove
            if(is_master()) print*, zk(k)
            qv(k) = qv_input(k)
  
          enddo
    end subroutine get_sounding

     real function interp_log( v_in, p_in, p_out, nzmax,nz_in  )
     implicit none
     integer, intent(in) ::  nz_in, nzmax
     real, intent(in) ::   v_in(nzmax), p_in(nzmax)
     real, intent(in) ::   p_out
    
     integer kp, k, im, ip, nz_out
     logical interp, increasing_z
     real    pres
     double precision :: w1, w2, v1, v2
     logical debug
     parameter ( debug = .false. )
    
     pres = p_out
    
     IF (pres > p_in(nz_in)) then
    !   if (is_master()) print*,"interp_log 1: p_in(nz_in), p_in(1) = " , p_in(nz_in), p_in(1)
       w2 = (log(p_in(nz_in))-log(pres))/(log(p_in(nz_in))-log(p_in(nz_in-1)))
       w1 = 1.-w2
       interp_log = v_in(nz_in)**w1 * v_in(nz_in-1)**w2
     ELSE IF (pres < p_in(1)) then ! extrapolate to lower pressure
       w2 = (log(p_in(2))-log(pres))/(log(p_in(2))-log(p_in(1)))
       w1 = 1.-w2
       v1 = v_in(1)
       v2 = v_in(2)
      ! interp_log =  w1*v_in(2) + w2*v_in(1) !
       interp_log = dble(v_in(2))**w1 * dble(v_in(1))**w2
    !  if (is_master()) print*,"interp_log 2: p_in(nz_in), p_in(1) = " , p_in(nz_in), p_in(1),pres,  &
    !              log(p_in(2))-log(pres), log(p_in(2))-log(p_in(1)), w2, w1, v_in(2) , v_in(1), &
    !              v_in(2)**w1, v_in(1)**w2, v2**w1 * v1**w2,  w1*v_in(2) + w2*v_in(1), interp_log
     ELSE
       interp = .false.
       kp = nz_in
       DO WHILE ( (interp .eqv. .false.) .and. (kp .ge. 2) )
         IF(   ((p_in(kp)   .ge. pres) .and.     &
                (p_in(kp-1) .le. pres))             )   THEN
           w2 = (log(pres)-log(p_in(kp)))/(log(p_in(kp-1))-log(p_in(kp)))
           w1 = 1.-w2
           interp_log = v_in(kp)**w1 * v_in(kp-1)**w2
           interp = .true.
    !    if (is_master()) print*,"interp_log 3: pres,p(kp),p(kp-1),w2,w1 = " , pres, p_in(kp), p_in(kp-1), &
    !               w2, w1, v_in(kp) , v_in(kp-1), interp_log,  w1*v_in(kp) + w2*v_in(kp-1)
    
         END IF
         kp = kp-1
       ENDDO
     ENDIF
    
     end function interp_log
    
     real function interp_lin( v_in, z_in, z_out, nzmax,nz_in  )
     implicit none
     integer, intent(in) ::  nz_in, nzmax
     real, intent(in) ::   v_in(nzmax), z_in(nzmax)
     real, intent(in) ::   z_out
    
     integer kp, k, im, ip, nz_out
     logical interp, increasing_z
     real    height, w1, w2
     logical debug
     parameter ( debug = .true. )
    
    height = z_out
    IF (height < z_in(nz_in)) then
          if(debug .and. is_master()) print*, ' point 1 in interp_lin ', height
          w2 = (z_in(nz_in)-height)/(z_in(nz_in)-z_in(nz_in-1))
          w1 = 1.-w2
          interp_lin = w1*v_in(nz_in) + w2*v_in(nz_in-1)
        ELSE IF (height > z_in(1)) then
          if(debug .and. is_master()) print*, ' point 2 in interp_lin ', height
          w2 = (z_in(2)-height)/(z_in(2)-z_in(1))
          w1 = 1.-w2
          interp_lin = w1*v_in(2) + w2*v_in(1)
        ELSE
          if(debug .and. is_master()) print*, ' point 3 in interp_lin ', height
          interp = .false.
          kp = nz_in
          height = z_out
          DO WHILE ( (interp .eqv. .false.) .and. (kp .ge. 2) )
            IF(   ((z_in(kp)   .le. height) .and.     &
                   (z_in(kp-1) .ge. height))             )   THEN
              w2 = (height-z_in(kp))/(z_in(kp-1)-z_in(kp))
              w1 = 1.-w2
              interp_lin = w1*v_in(kp) + w2*v_in(kp-1)
              interp = .true.
            END IF
            kp = kp-1
          ENDDO
        ENDIF
     end function interp_lin

     subroutine SuperCell_Sounding(km, ps, pk1, tp, qp)
       ! Morris Weisman & J. Klemp 2002 sounding
       ! Output sounding on pressure levels:
        integer, intent(in):: km
        real, intent(in):: ps     ! surface pressure (Pa)
        real, intent(in), dimension(km):: pk1
        real, intent(out), dimension(km):: tp, qp
       ! Local:
        integer, parameter:: ns = 401
        integer, parameter:: nx = 3
        real, dimension(ns):: zs, pt, qs, us, rh, pp, pk, dpk, dqdt
        real, parameter:: Tmin = 175.
        real, parameter:: p00 = 1.0e5
        real, parameter:: qst = 3.0e-6
        real, parameter:: qv0 = 1.4e-2
        real, parameter:: ztr = 12.E3
        real, parameter:: ttr = 213.
        real, parameter:: ptr = 343.    ! Tropopause potential temp.
        real, parameter:: pt0 = 300.    ! surface potential temperature
        real:: dz0, zvir, fac_z, pk0, temp1, p2
        integer:: k, n, kk
       
       ! tp=0.
       ! qp=0.
       
       !!#else
       
        zvir = rvgas/rdgas - 1.
        pk0 = p00**kappa
        pp(ns) = ps
        pk(ns) = ps**kappa
        if ( (is_master()) ) then
            write(*,*) 'Computing sounding for super-cell test'
        endif
       
        dz0 = 50.
        zs(ns) = 0.
        qs(:) = qst
        rh(:) = 0.25
       
        do k=ns-1, 1, -1
           zs(k) = zs(k+1) + dz0
        enddo
       
        do k=1,ns
       ! Potential temperature
           if ( zs(k) .gt. ztr ) then
       ! Stratosphere:
                pt(k) = ptr*exp(grav*(zs(k)-ztr)/(cp_air*ttr))
           else
       ! Troposphere:
                fac_z = (zs(k)/ztr)**1.25
                pt(k) = pt0 + (ptr-pt0)* fac_z
                rh(k) =  1. -     0.75 * fac_z
       ! First guess on q:
                qs(k) = qv0 - (qv0-qst)*fac_z
           endif
           pt(k) = pt(k) / pk0
        enddo
       
       !--------------------------------------
       ! Iterate nx times with virtual effect:
       !--------------------------------------
        do n=1, nx
           do k=1,ns-1
               temp1 = 0.5*(pt(k)*(1.+zvir*qs(k)) + pt(k+1)*(1.+zvir*qs(k+1)))
              dpk(k) = grav*(zs(k)-zs(k+1))/(cp_air*temp1)   ! DPK > 0
           enddo
       
           do k=ns-1,1,-1
              pk(k) = pk(k+1) - dpk(k)
           enddo
       
           do k=1, ns
              temp1 = pt(k)*pk(k)
       !      if ( (is_master()) ) write(*,*) k, temp1, rh(k)
              if ( pk(k) > 0. ) then
                   pp(k) = exp(log(pk(k))/kappa)
                   qs(k) = 380./pp(k)*exp(17.27*(temp1-273.)/(temp1-36.))
                   qs(k) = min( qv0, rh(k)*qs(k) )
                   if ( (is_master()) ) write(*,*) 0.01*pp(k), qs(k)
              else
                   if ( (is_master()) ) write(*,*) n, k, pk(k)
                   call mpp_error(FATAL, 'Super-Cell case: pk < 0')
              endif
           enddo
        enddo
       
       ! Interpolate to p levels using pk1: p**kappa
        do 555 k=1, km
           if ( pk1(k) .le. pk(1) ) then
                tp(k) = pt(1)*pk(1)/pk1(k)   ! isothermal above
                qp(k) = qst                  ! set to stratosphere value
           elseif ( pk1(k) .ge. pk(ns) ) then
                tp(k) = pt(ns)
                qp(k) = qs(ns)
           else
             do kk=1,ns-1
                if( (pk1(k).le.pk(kk+1)) .and. (pk1(k).ge.pk(kk)) ) then
                    fac_z = (pk1(k)-pk(kk))/(pk(kk+1)-pk(kk))
                    tp(k) = pt(kk) + (pt(kk+1)-pt(kk))*fac_z
                    qp(k) = qs(kk) + (qs(kk+1)-qs(kk))*fac_z
                    goto 555
                endif
             enddo
           endif
       555  continue
       
        do k=1,km
           tp(k) = tp(k)*pk1(k)    ! temperature
           tp(k) = max(Tmin, tp(k))
        enddo


     end subroutine SuperCell_Sounding

    end module fv_ideal_mod
