 module fv_surf_map_mod

      use fms_mod,           only: file_exist, check_nml_error,            &
                                   open_namelist_file, close_file, stdlog, &
                                   mpp_pe, mpp_root_pe, FATAL, error_mesg
      use mpp_mod,           only: get_unit, input_nml_file
      use mpp_domains_mod,   only: mpp_update_domains
      use constants_mod,     only: grav
#ifdef MARS_GCM
      use fms_mod,           only: read_data
      use fms_io_mod,        only: field_size
#endif

      use fv_grid_utils_mod, only: great_circle_dist, latlon2xyz, v_prod,  &
                                   sina_u, sina_v, g_sum, global_mx 
      use fv_mp_mod,         only: domain, ng, is,js,ie,je, isd,jsd,ied,jed, &
                                   gid, mp_stop, mp_reduce_min, mp_reduce_max
      use fv_timing_mod,     only: timing_on, timing_off

      implicit none

      real pi
      private
      real, allocatable:: sgh_g(:,:), oro_g(:,:), zs_g(:,:)
!-----------------------------------------------------------------------
! NAMELIST
!    Name, resolution, and format of XXmin USGS datafile
!      1min
!         nlon = 10800 * 2
!         nlat =  5400 * 2
!      2min
!         nlon = 10800
!         nlat =  5400
!      5min
!         nlon = 4320
!         nlat = 2160
!    surf_format:      netcdf (default)
!                      binary
      integer           ::  nlon = 10800
      integer           ::  nlat =  5400
#ifdef MARS_GCM
      character(len=128)::  surf_file = "INPUT/mars_topo.nc"
      character(len=6)  ::  surf_format = 'netcdf'
      character(len=80) :: field_name 
      integer           :: fld_dims(4)
      real, allocatable :: rtopo(:,:)
#else
      character(len=128)::  surf_file = "INPUT/topo5min.nc"
      character(len=6)  ::  surf_format = 'netcdf'
#endif
      namelist /surf_map_nml/ surf_file,surf_format,nlon,nlat
!
      public  sgh_g, oro_g, zs_g
      public  surfdrv, map_to_cubed_simple

      contains

      subroutine surfdrv(npx, npy, grid, agrid,   &
                         area, dx, dy, dxc, dyc, phis, master) 

      implicit         none
#include <netcdf.inc>
      integer, intent(in):: npx, npy
      logical master

    ! INPUT arrays
      real, intent(in)::area(is-ng:ie+ng, js-ng:je+ng)
      real, intent(in):: dx(is-ng:ie+ng, js-ng:je+ng+1)
      real, intent(in):: dy(is-ng:ie+ng+1, js-ng:je+ng)
      real, intent(in)::dxc(is-ng:ie+ng+1, js-ng:je+ng)
      real, intent(in)::dyc(is-ng:ie+ng, js-ng:je+ng+1)

      real, intent(in):: grid(is-ng:ie+ng+1, js-ng:je+ng+1,2)
      real, intent(in):: agrid(is-ng:ie+ng, js-ng:je+ng,2)

    ! OUTPUT arrays
      real, intent(out):: phis(is-ng:ie+ng, js-ng:je+ng)
! Local:
      real :: z2(is:ie, js:je)
! Position of edges of the box containing the original data point:

      integer          londim
      integer          latdim

      real dx1, dx2, dy1, dy2

      character(len=80) :: topoflnm
      real(kind=4) :: fmin, fmax
      real(kind=4), allocatable :: ftopo(:,:), htopo(:,:)
      real, allocatable :: lon1(:),  lat1(:)
      integer i, j, n
      integer ncid, lonid, latid, ftopoid, htopoid
      integer status
      logical check_orig
      real da_min, da_max, cd2, cd4, zmean, z2mean
      integer fid

! Output the original 10 min NAVY data in grads readable format
      data             check_orig /.false./

      allocate ( oro_g(isd:ied, jsd:jed) )
      allocate ( sgh_g(isd:ied, jsd:jed) )
      allocate (  zs_g(is:ie, js:je) )


      call read_namelist

#ifdef MARS_GCM
      if (surf_format == "binary")  &
          call error_mesg ( 'surfdrv', ' binary input not allowed for Mars Model', FATAL )

      if (file_exist(surf_file)) then
         field_name = 'topo'
 
!         call field_size( trim(surf_file), 'lat', fld_dims )
!         call field_size( trim(surf_file), 'lon', fld_dims )
         call field_size( trim(surf_file), trim(field_name), fld_dims )

         nlon= fld_dims(1);  nlat= fld_dims(2);

         if(master) write(*,*) 'Mars Terrain dataset dims=', nlon, nlat

         allocate ( htopo(nlon,nlat) )
         allocate ( rtopo(nlon,nlat) )
         allocate ( ftopo(nlon,nlat) )

         call read_data( trim(surf_file), trim(field_name), rtopo, no_domain=.true. )

!   This is needed because htopo is declared as real*4
         htopo= rtopo
         ftopo = 1.0              ! all land points

         if ( master ) then
            write(6,*) 'Check Hi-res Mars data ..'
            fmax =  vmax(htopo,fmin,nlon,nlat,1)
            write(6,*) 'hmax=', fmax
            write(6,*) 'hmin=', fmin
         endif
      else
         call error_mesg ( 'surfdrv', 'mars_topo not found in INPUT', FATAL )
      endif
#else

      if (file_exist(surf_file)) then
!
! surface file in NetCDF format
!
        if (surf_format == "netcdf") then

          allocate ( ftopo(nlon,nlat) )
          allocate ( htopo(nlon,nlat) )

          if ( master ) write(*,*) 'Opening USGS datset file:', surf_file, surf_format, nlon, nlat
  
          status = nf_open (surf_file, NF_NOWRITE, ncid)
          if (status .ne. NF_NOERR) call handle_err(status)
  
          status = nf_inq_dimid (ncid, 'lon', lonid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_inq_dimlen (ncid, lonid, londim)
          if (status .ne. NF_NOERR) call handle_err(status)

          status = nf_inq_dimid (ncid, 'lat', latid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_inq_dimlen (ncid, latid, latdim)
          if (status .ne. NF_NOERR) call handle_err(status)

          status = nf_inq_varid (ncid, 'ftopo', ftopoid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_get_var_real (ncid, ftopoid, ftopo)
          if (status .ne. NF_NOERR) call handle_err(status)

          status = nf_inq_varid (ncid, 'htopo', htopoid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_get_var_real (ncid, htopoid, htopo)
          if (status .ne. NF_NOERR) call handle_err(status)

          status = nf_close (ncid)
          if (status .ne. NF_NOERR) call handle_err(status)
!
! ... Set check_orig=.true. to output original 10-minute
!         real(kind=4) ::  data (GrADS format)
!
          if (check_orig) then
            topoflnm = 'topo.bin'
            fid = get_unit()
            open (unit=fid, file=topoflnm, form='unformatted', &
                status='unknown', access='direct', recl=nlon*nlat*4)
            write (fid, rec=1) ftopo
            write (fid, rec=2) htopo
            close (unit=fid)
          endif
  
!
! surface file in binary format
!
        elseif (surf_format == "binary") then 
!        nlon = 10800
!        nlat =  5400
!        surf_file    = '/work/sjl/topo/topo2min.bin'
  
          if ( master ) write(*,*) 'Opening USGS datset file:', surf_file, surf_format, nlon, nlat

          fid = get_unit()
          open (unit=fid, file=surf_file, form='unformatted', &
                    status='unknown', access='direct', recl=nlon*nlat*4)

          allocate ( ftopo(nlon,nlat) )
          allocate ( htopo(nlon,nlat) )
          read (fid, rec=1) ftopo
          read (fid, rec=2) htopo
          close (unit=fid)
        endif

      else
        if(master) write(*,*) 'USGS dataset = ', surf_file, surf_format
        call error_mesg ( 'surfdrv',  &
            'missing input file', FATAL )
      endif
#endif MARS_GCM

      allocate ( lat1(nlat+1) )
      allocate ( lon1(nlon+1) )

      pi = 4.0 * datan(1.0d0)

      dx1 = 2.*pi/real(nlon)
      dy1 = pi/real(nlat)

      do i=1,nlon+1
         lon1(i) = dx1 * real(i-1)    ! between 0 2pi
      enddo

         lat1(1) = - 0.5*pi
         lat1(nlat+1) =  0.5*pi
      do j=2,nlat
         lat1(j) = -0.5*pi + dy1*(j-1)
      enddo

      if ( master ) then
           write(6,*) 'check original data ..'
           fmax =  vmax(htopo,fmin,nlon,nlat,1)
           write(6,*) 'hmax=', fmax
           write(6,*) 'hmin=', fmin
      endif

!-------------------------------------
! Compute raw phis and oro
!-------------------------------------

                                                      call timing_on('map_to_cubed')
      call map_to_cubed_raw(nlon, nlat, lat1, lon1, htopo, ftopo, grid, agrid,  &
                            phis, oro_g, sgh_g, master, npx, npy)

      if(master) write(*,*) 'map_to_cubed_raw: done'
!     write(*,*) gid, 'map_to_cubed_raw: done'
                                                      call timing_off('map_to_cubed')
      deallocate ( htopo )
      deallocate ( ftopo )
      deallocate ( lon1 )
      deallocate ( lat1 )

#ifndef MARS_GCM
      call remove_ice_sheets (agrid(isd,jsd,1), agrid(isd,jsd,2), oro_g )
#endif

      call global_mx(oro_g, ng, da_min, da_max)
      if ( master ) write(*,*) 'ORO min=', da_min, ' Max=', da_max

      do j=js,je
         do i=is,ie
            zs_g(i,j) = phis(i,j)
            z2(i,j) = phis(i,j)**2
         end do
      end do
!--------
! Filter:
!--------
      call global_mx(phis, ng, da_min, da_max)
      zmean = g_sum(zs_g, is, ie, js, je, ng, area, 1)
      z2mean = g_sum(z2 , is, ie, js, je, ng, area, 1)

      if ( master ) then
           write(*,*) 'Before filter ZS min=', da_min, ' Max=', da_max,' Mean=',zmean
           write(*,*) '*** Mean variance *** =', z2mean
      endif

      call global_mx(area, ng, da_min, da_max)

                                                    call timing_on('Terrain_filter')
! Del-2:
      cd2 = 0.20*da_min
      if ( npx>=721 ) then
           if ( npx<=1001 ) then
                call del2_cubed_sphere(npx, npy, phis, area, dx, dy, dxc, dyc, 1, cd2)
           else
                call del2_cubed_sphere(npx, npy, phis, area, dx, dy, dxc, dyc, 2, cd2)
           endif
      endif

! MFCT Del-4:
      if ( npx<=91 ) then
         cd4 = 0.20 * da_min
         call del4_cubed_sphere(npx, npy, phis, area, dx, dy, dxc, dyc, 1, cd4)
      elseif( npx<=181 ) then
         cd4 = 0.20 * da_min
         call del4_cubed_sphere(npx, npy, phis, area, dx, dy, dxc, dyc, 2, cd4)
      elseif( npx<=361 ) then
         cd4 = 0.20 * da_min
         call del4_cubed_sphere(npx, npy, phis, area, dx, dy, dxc, dyc, 4, cd4)
      elseif( npx<=721 ) then
         cd4 = 0.20 * da_min
         call del4_cubed_sphere(npx, npy, phis, area, dx, dy, dxc, dyc, 6, cd4)
      else
         cd4 = 0.20 * da_min
         call del4_cubed_sphere(npx, npy, phis, area, dx, dy, dxc, dyc, 8, cd4)
      endif

      do j=js,je
         do i=is,ie
            z2(i,j) = phis(i,j)**2
         end do
      end do

      call global_mx(phis, ng, da_min, da_max)
      zmean  = g_sum(phis(is:ie,js:je), is, ie, js, je, ng, area, 1)
      z2mean = g_sum(z2,                is, ie, js, je, ng, area, 1)

      if ( master ) then
           write(*,*) 'After filter Phis min=', da_min, ' Max=', da_max, 'Mean=', zmean
           write(*,*) '*** Mean variance *** =', z2mean
      endif

      do j=js,je
         do i=is,ie
            phis(i,j) =  grav * phis(i,j)
            if ( sgh_g(i,j) <= 0. ) then
                 sgh_g(i,j) = 0.
            else
                 sgh_g(i,j) = sqrt(sgh_g(i,j))
            endif
#ifdef SET_FLAG
            if ( oro_g(i,j) > .5 ) then
                oro_g(i,j) = 1.
            else
                oro_g(i,j) = 0.
            endif
#endif
         end do
      end do

      call global_mx(sgh_g, ng, da_min, da_max)
      if ( master ) write(*,*) 'Before filter SGH min=', da_min, ' Max=', da_max


!-----------------------------------------------
! Filter the standard deviation of mean terrain:
!-----------------------------------------------
      call global_mx(area, ng, da_min, da_max)
      call del4_cubed_sphere(npx, npy, sgh_g, area, dx, dy, dxc, dyc, 1, cd4)
      call global_mx(sgh_g, ng, da_min, da_max)
      if ( master ) write(*,*) 'After filter SGH min=', da_min, ' Max=', da_max
      do j=js,je
         do i=is,ie
            sgh_g(i,j) = max(0., sgh_g(i,j))
         enddo
      enddo
                                                    call timing_off('Terrain_filter')


 end subroutine surfdrv



 subroutine del2_cubed_sphere(npx, npy, q, area, dx, dy, dxc, dyc, nmax, cd)
      integer, intent(in):: npx, npy
      integer, intent(in):: nmax
      real, intent(in):: cd
    ! INPUT arrays
      real, intent(in)::area(isd:ied,  jsd:jed)
      real, intent(in)::  dx(isd:ied,  jsd:jed+1)
      real, intent(in)::  dy(isd:ied+1,jsd:jed)
      real, intent(in):: dxc(isd:ied+1,jsd:jed)
      real, intent(in):: dyc(isd:ied,  jsd:jed+1)
    ! OUTPUT arrays
      real, intent(inout):: q(is-ng:ie+ng, js-ng:je+ng)
! Local:
      real ddx(is:ie+1,js:je), ddy(is:ie,js:je+1)
      integer i,j,n

      call mpp_update_domains(q,domain,whalo=ng,ehalo=ng,shalo=ng,nhalo=ng)

! First step: average the corners:
      if ( is==1 .and. js==1 ) then
           q(1,1) = (q(1,1)+q(0,1)+q(1,0)) / 3.
           q(0,1) =  q(1,1)
           q(1,0) =  q(1,1)
      endif
      if ( (ie+1)==npx .and. js==1 ) then
           q(ie, 1) = (q(ie,1)+q(npx,1)+q(ie,0)) / 3.
           q(npx,1) =  q(ie,1)
           q(ie, 0) =  q(ie,1)
      endif
      if ( (ie+1)==npx .and. (je+1)==npy ) then
           q(ie, je) = (q(ie,je)+q(npx,je)+q(ie,npy)) / 3.
           q(npx,je) =  q(ie,je)
           q(ie,npy) =  q(ie,je)
      endif
      if ( is==1 .and. (je+1)==npy ) then
           q(1, je) = (q(1,je)+q(0,je)+q(1,npy)) / 3.
           q(0, je) =  q(1,je)
           q(1,npy) =  q(1,je)
      endif


      do n=1,nmax

         if( n>1 ) call mpp_update_domains(q,domain,whalo=ng,ehalo=ng,shalo=ng,nhalo=ng)

         do j=js,je
            do i=is,ie+1
               ddx(i,j) = dy(i,j)*sina_u(i,j)*(q(i-1,j)-q(i,j))/dxc(i,j)
            enddo
         enddo

         do j=js,je+1
            do i=is,ie
               ddy(i,j) = dx(i,j)*sina_v(i,j)*(q(i,j-1)-q(i,j))/dyc(i,j)
            enddo
         enddo

         do j=js,je
            do i=is,ie
               q(i,j) = q(i,j) + cd/area(i,j)*(ddx(i,j)-ddx(i+1,j)+ddy(i,j)-ddy(i,j+1))
            enddo
         enddo
       
      enddo

 end subroutine del2_cubed_sphere


 subroutine del4_cubed_sphere(npx, npy, q, area, dx, dy, dxc, dyc, nmax, cd)
      real, parameter:: esl = 1.E-20
      integer, intent(in):: npx, npy, nmax
      real, intent(in)::area(isd:ied,  jsd:jed)
      real, intent(in)::  dx(isd:ied,  jsd:jed+1)
      real, intent(in)::  dy(isd:ied+1,jsd:jed)
      real, intent(in):: dxc(isd:ied+1,jsd:jed)
      real, intent(in):: dyc(isd:ied,  jsd:jed+1)
      real, intent(in):: cd
      real, intent(inout):: q(is-ng:ie+ng, js-ng:je+ng)
! diffusive fluxes: 
      real :: fx2(is:ie+1,js:je), fy2(is:ie,js:je+1)
      real :: fx4(is:ie+1,js:je), fy4(is:ie,js:je+1)
      real   d2(isd:ied,jsd:jed)
      real  win(isd:ied,jsd:jed)
      real  wou(isd:ied,jsd:jed)

      real qlow(is:ie,js:je)
      real qmin(is:ie,js:je)
      real qmax(is:ie,js:je)
      integer i,j, n

  do n=1,nmax
      call mpp_update_domains(q,domain)

! First step: average the corners:
      if ( is==1 .and. js==1 ) then
           q(1,1) = (q(1,1)+q(0,1)+q(1,0)) / 3.
           q(0,1) =  q(1,1)
           q(1,0) =  q(1,1)
      endif
      if ( (ie+1)==npx .and. js==1 ) then
           q(ie, 1) = (q(ie,1)+q(npx,1)+q(ie,0)) / 3.
           q(npx,1) =  q(ie,1)
           q(ie, 0) =  q(ie,1)
      endif
      if ( (ie+1)==npx .and. (je+1)==npy ) then
           q(ie, je) = (q(ie,je)+q(npx,je)+q(ie,npy)) / 3.
           q(npx,je) =  q(ie,je)
           q(ie,npy) =  q(ie,je)
      endif
      if ( is==1 .and. (je+1)==npy ) then
           q(1, je) = (q(1,je)+q(0,je)+q(1,npy)) / 3.
           q(0, je) =  q(1,je)
           q(1,npy) =  q(1,je)
      endif

     do j=js,je
        do i=is,ie
           qmin(i,j) = min(q(i,j-1), q(i-1,j), q(i,j), q(i+1,j), q(i,j+1))
           qmax(i,j) = max(q(i,j-1), q(i-1,j), q(i,j), q(i+1,j), q(i,j+1))
        enddo
     enddo

!--------------
! Compute del-2
!--------------
!     call copy_corners(q, npx, npy, 1)
      do j=js,je
         do i=is,ie+1
            fx2(i,j) = cd*dy(i,j)*sina_u(i,j)*(q(i-1,j)-q(i,j))/dxc(i,j)
         enddo
      enddo

!     call copy_corners(q, npx, npy, 2)
      do j=js,je+1
         do i=is,ie
            fy2(i,j) = cd*dx(i,j)*sina_v(i,j)*(q(i,j-1)-q(i,j))/dyc(i,j)
         enddo
      enddo

      do j=js,je
         do i=is,ie
            d2(i,j) = (fx2(i,j)-fx2(i+1,j)+fy2(i,j)-fy2(i,j+1)) / area(i,j)
! Low order monotonic solution
            qlow(i,j) = q(i,j) + d2(i,j)
            d2(i,j) = cd * d2(i,j)
         enddo
      enddo

      call mpp_update_domains(d2,domain)

!---------------------
! Compute del4 fluxes:
!---------------------
!     call copy_corners(d2, npx, npy, 1)
      do j=js,je
         do i=is,ie+1
            fx4(i,j) = dy(i,j)*sina_u(i,j)*(d2(i,j)-d2(i-1,j))/dxc(i,j)-fx2(i,j)
         enddo
      enddo

!     call copy_corners(d2, npx, npy, 2)
      do j=js,je+1
         do i=is,ie
            fy4(i,j) = dx(i,j)*sina_v(i,j)*(d2(i,j)-d2(i,j-1))/dyc(i,j)-fy2(i,j)
         enddo
      enddo

!----------------
! Flux limitting:
!----------------
#ifndef NO_MFCT_FILTER
      do j=js,je
         do i=is,ie
            win(i,j) = max(0.,fx4(i,  j)) - min(0.,fx4(i+1,j)) +   &
                       max(0.,fy4(i,  j)) - min(0.,fy4(i,j+1)) + esl
            wou(i,j) = max(0.,fx4(i+1,j)) - min(0.,fx4(i,  j)) +   &
                       max(0.,fy4(i,j+1)) - min(0.,fy4(i,  j)) + esl
            win(i,j) = max(0., qmax(i,j) - qlow(i,j)) / win(i,j)*area(i,j)
            wou(i,j) = max(0., qlow(i,j) - qmin(i,j)) / wou(i,j)*area(i,j)
         enddo
      enddo

      call mpp_update_domains(win,domain, complete=.false.)
      call mpp_update_domains(wou,domain, complete=.true.)

      do j=js,je
         do i=is,ie+1
            if ( fx4(i,j) > 0. ) then
                 fx4(i,j) = min(1., wou(i-1,j), win(i,j)) * fx4(i,j) 
            else
                 fx4(i,j) = min(1., win(i-1,j), wou(i,j)) * fx4(i,j) 
            endif
         enddo
      enddo
      do j=js,je+1
         do i=is,ie
            if ( fy4(i,j) > 0. ) then
                 fy4(i,j) = min(1., wou(i,j-1), win(i,j)) * fy4(i,j) 
            else
                 fy4(i,j) = min(1., win(i,j-1), wou(i,j)) * fy4(i,j) 
            endif
         enddo
      enddo
#endif

! Update:
      do j=js,je
         do i=is,ie
            q(i,j) = qlow(i,j) + (fx4(i,j)-fx4(i+1,j)+fy4(i,j)-fy4(i,j+1))/area(i,j)
         enddo
      enddo

  enddo    ! end n-loop

 end subroutine del4_cubed_sphere




 subroutine map_to_cubed_raw(im, jm, lat1, lon1, q1, f1,  grid, agrid,  &
                                  q2, f2, h2, master, npx, npy)

! Input
      integer, intent(in):: im, jm        ! original dimensions
      integer, intent(in):: npx, npy
      real, intent(in):: lat1(jm+1)       ! original southern edge of the cell [-pi/2:pi/2]
      real, intent(in):: lon1(im+1)       ! original western edge of the cell [0:2*pi]
      real(kind=4), intent(in):: q1(im,jm)      ! original data at center of the cell
      real(kind=4), intent(in):: f1(im,jm)      !

      real, intent(in)::  grid(isd:ied+1, jsd:jed+1,2)
      real, intent(in):: agrid(isd:ied,   jsd:jed,  2)
      logical, intent(in):: master
! Output
      real, intent(out):: q2(isd:ied,jsd:jed) ! Mapped data at the target resolution
      real, intent(out):: f2(isd:ied,jsd:jed) ! oro
      real, intent(out):: h2(isd:ied,jsd:jed) ! variances of terrain

! Local
      real(kind=4), allocatable:: qt(:,:), ft(:,:), lon_g(:)
      real lat_g(jm)
      real pc(3), p2(2), pp(3), grid3(3, is:ie+1, js:je+1)
      integer i,j, np
      integer igh
      integer ii, jj, i1, i2, j1, j2
      integer ifirst, ilast
      real qsum, fsum, hsum, lon_w, lon_e, lat_s, lat_n, r2d
      real delg, dlat
!     integer, parameter:: lat_crit = 15             ! 15 * (1/30) = 0.5 deg
      integer:: lat_crit
      integer, parameter:: ig = 2
      real q1_np, q1_sp, f1_np, f1_sp, h1_sp, h1_np, pi5, deg0
      logical inside

      pi5 = 0.5 * pi
      r2d = 180./pi

!     lat_crit = jm / min(360, 4*(npx-1))    ! 0.5  (deg) or larger
      lat_crit = jm / min(720, 8*(npx-1))    ! 0.25 (deg) or larger

      dlat = 180./real(jm)

      igh = im/4 + 1

      if (master) write(*,*) 'Terrain dataset im=', im, 'jm=', jm
      if (master) write(*,*) 'igh (terrain ghosting)=', igh

      allocate (    qt(-igh:im+igh,jm) )
      allocate (    ft(-igh:im+igh,jm) )
      allocate ( lon_g(-igh:im+igh   ) )

! Ghost the input coordinates:
      do i=1,im
         lon_g(i) = 0.5*(lon1(i)+lon1(i+1))
      enddo

      do i=-igh,0
         lon_g(i) = lon_g(i+im)
      enddo
      do i=im+1,im+igh
         lon_g(i) = lon_g(i-im)
      enddo

      do j=1,jm
         lat_g(j) = 0.5*(lat1(j)+lat1(j+1))
      enddo

!     if ( 2*(im/2) /= im ) then
!          write(*,*) 'Warning: Terrain datset must have an even nlon dimension'
!     endif

! Ghost Data
      do j=1,jm
         do i=1,im
            qt(i,j) = q1(i,j)
            ft(i,j) = f1(i,j)
         enddo
         do i=-igh,0
            qt(i,j) = qt(i+im,j)
            ft(i,j) = ft(i+im,j)
         enddo
         do i=im+1,im+igh
            qt(i,j) = qt(i-im,j)
            ft(i,j) = ft(i-im,j)
         enddo
      enddo

      do j=js,je+1
         do i=is,ie+1
            call latlon2xyz(grid(i,j,1:2), grid3(1,i,j))
         enddo
      enddo

! Compute values very close to the poles:
!----
! SP:
!----
     qsum = 0.
     fsum = 0.
     hsum = 0.
     np   = 0
     do j=1,lat_crit
        do i=1,im
           np = np + 1
           qsum = qsum + q1(i,j)
           fsum = fsum + f1(i,j)
        enddo
     enddo
     q1_sp = qsum / real(np)
     f1_sp = fsum / real(np)

     hsum = 0.
     do j=1,lat_crit
        do i=1,im
           hsum = hsum + (q1_sp-q1(i,j))**2
        enddo
     enddo
     h1_sp = hsum / real(np)

     if(master) write(*,*) 'SP:', q1_sp, f1_sp, sqrt(h1_sp)
!----
! NP:
!----
     qsum = 0.
     fsum = 0.
     hsum = 0.
     np   = 0
     do j=jm-lat_crit+1,jm
        do i=1,im
           np = np + 1
           qsum = qsum + q1(i,j)
           fsum = fsum + f1(i,j)
        enddo
     enddo
     q1_np = qsum / real(np)
     f1_np = fsum / real(np)

     hsum = 0.
     do j=jm-lat_crit+1,jm
        do i=1,im
           hsum = hsum + (q1_np-q1(i,j))**2
        enddo
     enddo
     h1_np = hsum / real(np)

     if(master) write(*,*) 'NP:', q1_np, f1_np, sqrt(h1_np)
     if(master) write(*,*) 'surf_map: Search started ....'

      do 4444 j=js,je
         do 4444 i=is,ie
 
            lat_s = min( grid(i,j,2), grid(i+1,j,2), grid(i,j+1,2), grid(i+1,j+1,2) )
            lat_n = max( grid(i,j,2), grid(i+1,j,2), grid(i,j+1,2), grid(i+1,j+1,2) )

            if ( r2d*lat_n < (lat_crit*dlat - 90.) ) then
                 q2(i,j) = q1_sp
                 f2(i,j) = f1_sp
                 h2(i,j) = h1_sp
                 go to 4444
            elseif ( r2d*lat_s > (90. - lat_crit*dlat) ) then
                 q2(i,j) = q1_np
                 f2(i,j) = f1_np
                 h2(i,j) = h1_np
                 go to 4444
            endif

            j1 = nint( (pi5+lat_s)/(pi/real(jm)) ) - ig
            j2 = nint( (pi5+lat_n)/(pi/real(jm)) ) + ig
            j1 = max(1,  j1)
            j2 = min(jm, j2)

            lon_w = min( grid(i,j,1), grid(i+1,j,1), grid(i,j+1,1), grid(i+1,j+1,1) ) 
            lon_e = max( grid(i,j,1), grid(i+1,j,1), grid(i,j+1,1), grid(i+1,j+1,1) )
            if ( (lon_e - lon_w) > pi ) then
                 i1 = nint( (lon_e-2.*pi)/(2.*pi/real(im)) )
                 i2 = nint(  lon_w       /(2.*pi/real(im)) )
            else
                 i1 = nint( lon_w / (2.*pi/real(im)) )
                 i2 = nint( lon_e / (2.*pi/real(im)) )
            endif

            i1 = max(  -igh, i1 - ig)
            i2 = min(im+igh, i2 + ig)

              np = 0
            qsum = 0.
            fsum = 0.
            hsum = 0.
            do jj=j1,j2
               p2(2) = lat_g(jj)
               do ii=i1,i2
                  p2(1) = lon_g(ii)
                  call latlon2xyz(p2, pp)
                  inside=inside_p4(grid3(1,i,j), grid3(1,i+1,j), grid3(1,i+1,j+1), grid3(1,i,j+1), pp)
                  if ( inside ) then
                      np = np + 1
                      qsum = qsum + qt(ii,jj)
                      fsum = fsum + ft(ii,jj)
                      hsum = hsum + qt(ii,jj)**2
                  endif
               enddo
            enddo

            if ( np > 0 ) then
                 q2(i,j) = qsum / real(np)
                 f2(i,j) = fsum / real(np)
                 h2(i,j) = hsum / real(np) - q2(i,j)**2
            else
                 write(*,*) 'Surf_map failed for GID=', gid, '(lon,lat)=', agrid(i,j,1)*r2d,agrid(i,j,2)*r2d
                 stop
!                call mp_stop   ! does not really stop !!!
            endif

4444  continue

      deallocate (   qt )
      deallocate (   ft )
      deallocate (lon_g )

 end subroutine map_to_cubed_raw



 logical function inside_p4(p1, p2, p3, p4, pp)
!
!            4----------3
!           /          /
!          /    pp    /
!         /          /
!        1----------2
!
! A * B = |A| |B| cos(angle)

      real, intent(in):: p1(3), p2(3), p3(3), p4(3)
      real, intent(in):: pp(3)
! Local:
      real v1(3), v2(3), vp(3)
      real a1, a2, aa, s1, s2, ss
      integer k

! S-W:
      do k=1,3
         v1(k) = p2(k) - p1(k) 
         v2(k) = p4(k) - p1(k) 
         vp(k) = pp(k) - p1(k) 
      enddo
      s1 = sqrt( v1(1)**2 + v1(2)**2 + v1(3)**2 )
      s2 = sqrt( v2(1)**2 + v2(2)**2 + v2(3)**2 )
      ss = sqrt( vp(1)**2 + vp(2)**2 + vp(3)**2 )

! Compute cos(angle):
      aa = v_prod(v1, v2) / (s1*s2)
      a1 = v_prod(v1, vp) / (s1*ss)
      a2 = v_prod(v2, vp) / (s2*ss)

      if ( a1<aa  .or.  a2<aa ) then
           inside_p4 = .false.
           return
      endif

! N-E:
      do k=1,3
         v1(k) = p2(k) - p3(k) 
         v2(k) = p4(k) - p3(k) 
         vp(k) = pp(k) - p3(k) 
      enddo
      s1 = sqrt( v1(1)**2 + v1(2)**2 + v1(3)**2 )
      s2 = sqrt( v2(1)**2 + v2(2)**2 + v2(3)**2 )
      ss = sqrt( vp(1)**2 + vp(2)**2 + vp(3)**2 )

! Compute cos(angle):
      aa = v_prod(v1, v2) / (s1*s2)
      a1 = v_prod(v1, vp) / (s1*ss)
      a2 = v_prod(v2, vp) / (s2*ss)

      if ( a1<aa  .or.  a2<aa ) then
           inside_p4 = .false.
      else
           inside_p4 = .true.
      endif

 end function inside_p4



 subroutine handle_err(status)
#include <netcdf.inc>
      integer          status

      if (status .ne. nf_noerr) then
        print *, nf_strerror(status)
        stop 'Stopped'
      endif

 end subroutine  handle_err


 real function vmax(a,pmin,m,n,z)
      integer m,n,z, i,j,k
      real(kind=4) :: pmin, pmax
      real(kind=4) :: a(m,n,z)

      pmax = a(1,1,1)
      pmin = a(1,1,1)

      do k=1,z
      do j=1,n
      do i=1,m
         pmax = max(pmax,a(i,j,k))
         pmin = min(pmin,a(i,j,k))
      enddo
      enddo
      enddo   

      vmax = pmax
 end function vmax

         
 subroutine remove_ice_sheets (lon, lat, lfrac )
!---------------------------------
! Bruce Wyman's fix for Antarctic
!--------------------------------- 
      real, intent(in)    :: lon(isd:ied,jsd:jed), lat(isd:ied,jsd:jed)
      real, intent(inout) :: lfrac(isd:ied,jsd:jed)
        
! lon   = longitude in radians
! lat   = latitude in radians
! lfrac = land-sea mask (land=1, sea=0)
            
      integer :: i, j
      real :: dtr, phs, phn
            
      dtr = acos(0.)/90.
      phs = -83.9999*dtr                                  
!     phn = -78.9999*dtr
      phn = -76.4*dtr
            
      do j = js, je
         do i = is, ie
         if ( lat(i,j) < phn ) then
                              ! replace all below this latitude
         if ( lat(i,j) < phs ) then
              lfrac(i,j) = 1.0
              cycle
         endif
                              ! replace between 270 and 360 deg
         if ( sin(lon(i,j)) < 0. .and. cos(lon(i,j)) > 0.) then
              lfrac(i,j) = 1.0
              cycle 
         endif
         endif
         enddo
      enddo
 end subroutine remove_ice_sheets

    subroutine map_to_cubed_simple(im, jm, lat1, lon1, q1, grid, agrid, q2, npx, npy)


! Input
      integer, intent(in):: im,jm         ! original dimensions
      integer, intent(in):: npx, npy
!rjw      logical, intent(in):: master
      real, intent(in):: lat1(jm+1)       ! original southern edge of the cell [-pi/2:pi/2]
      real, intent(in):: lon1(im+1)       ! original western edge of the cell [0:2*pi]
      real(kind=4), intent(in):: q1(im,jm)        ! original data at center of the cell
!rjw      real(kind=4), intent(in):: f1(im,jm)        !

      real, intent(in)::  grid(is-ng:ie+ng+1, js-ng:je+ng+1,2)
      real, intent(in):: agrid(is-ng:ie+ng,   js-ng:je+ng,  2)

! Output
      real, intent(out):: q2(is-ng:ie+ng, js-ng:je+ng) ! Mapped data at the target resolution
!rjw      real, intent(out):: f2(isd:ied,jsd:jed) ! oro
!rjw      real, intent(out):: h2(isd:ied,jsd:jed) ! variances of terrain

! Local
      real(kind=4)  qt(-im/32:im+im/32,jm)    ! ghosted east-west
!rjw      real(kind=4)  ft(-im/32:im+im/32,jm)    ! 
      real lon_g(-im/32:im+im/32)
      real lat_g(jm)

      real pc(3), p2(2), pp(3), grid3(3,is-ng:ie+ng+1, js-ng:je+ng+1)
      integer i,j, np
      integer ii, jj, i1, i2, j1, j2
      integer ifirst, ilast
      real ddeg, latitude, qsum, fsum, hsum, lon_w, lon_e, lat_s, lat_n, r2d
      real delg

      pi = 4.0 * datan(1.0d0)

      r2d = 180./pi
      ddeg = 2.*pi/real(4*npx)

! Ghost the input coordinates:
      do i=1,im
         lon_g(i) = 0.5*(lon1(i)+lon1(i+1))
      enddo

      do i=-im/32,0
         lon_g(i) = lon_g(i+im)
      enddo
      do i=im+1,im+im/32
         lon_g(i) = lon_g(i-im)
      enddo

      do j=1,jm
         lat_g(j) = 0.5*(lat1(j)+lat1(j+1))
      enddo

      if ( 2*(im/2) /= im ) then
           write(*,*) 'Warning: Terrain datset must have an even nlon dimension'
      endif
! Ghost Data
      do j=1,jm
         do i=1,im
            qt(i,j) = q1(i,j)
!rjw            ft(i,j) = f1(i,j)
         enddo
         do i=-im/32,0
            qt(i,j) = qt(i+im,j)
!rjw            ft(i,j) = ft(i+im,j)
         enddo
         do i=im+1,im+im/32
            qt(i,j) = qt(i-im,j)
!rjw            ft(i,j) = ft(i-im,j)
         enddo
      enddo
      
      do j=js,je+1
         do i=is,ie+1
            call latlon2xyz(grid(i,j,1:2), grid3(1,i,j))
         enddo
      enddo

!rjw     if(master) write(*,*) 'surf_map: Search started ....'
! Mapping:
      do j=js,je
         do i=is,ie
! Determine the approximate local loop bounds (slightly larger than needed)
            lon_w = min( grid(i,j,1), grid(i+1,j,1), grid(i,j+1,1), grid(i+1,j+1,1) ) - ddeg
            lon_e = max( grid(i,j,1), grid(i+1,j,1), grid(i,j+1,1), grid(i+1,j+1,1) ) + ddeg
            if ( (lon_e - lon_w) > pi ) then
                 delg = max( abs(lon_e-2.*pi), abs(lon_w) ) + ddeg
                 i1 = -delg / (2.*pi/real(im)) - 1
                 i2 = -i1 + 1
            else 
                 i1 = lon_w / (2.*pi/real(im)) - 1
                 i2 = lon_e / (2.*pi/real(im)) + 2
            endif
            i1 = max(-im/32, i1)
            i2 = min(im+im/32, i2)
!           
            lat_s = min( grid(i,j,2), grid(i+1,j,2), grid(i,j+1,2), grid(i+1,j+1,2) ) - ddeg
            lat_n = max( grid(i,j,2), grid(i+1,j,2), grid(i,j+1,2), grid(i+1,j+1,2) ) + ddeg
            j1 = (0.5*pi + lat_s) / (pi/real(jm)) - 1
            j2 = (0.5*pi + lat_n) / (pi/real(jm)) + 2
              
              np = 0
            qsum = 0.
!rjw            fsum = 0.
!rjw            hsum = 0.
!           call latlon2xyz(agrid(i,j,1:2), pc)

!rjw             print *, 'Interior loop:  ',  i, j, i1, i2, j1, j2,  grid(i,j,1:2)*r2d,  agrid(i,j,1:2)*r2d

            do jj=max(1,j1),min(jm,j2)
                  p2(2) = lat_g(jj)
                  latitude =  p2(2)*r2d
               if ( abs(latitude) > 80.  ) then
                  ifirst = 1; ilast = im
               else
                  ifirst = i1; ilast = i2
               endif

               do ii=ifirst, ilast
                  p2(1) = lon_g(ii)
                  call latlon2xyz(p2, pp)
                  if (inside_p4(grid3(1,i,j), grid3(1,i+1,j), grid3(1,i+1,j+1), grid3(1,i,j+1), pp)) then
                       np = np + 1
                       qsum = qsum + qt(ii,jj)
!rjw                       fsum = fsum + ft(ii,jj)
!rjw                       hsum = hsum + qt(ii,jj)**2
                  endif

               enddo
            enddo
! Compute weighted average:
            if ( np > 0 ) then
                 q2(i,j) = qsum / real(np)
!rjw                 f2(i,j) = fsum / real(np)
!rjw                 h2(i,j) = hsum / real(np) - q2(i,j)**2
            else                    ! the subdomain could be totally flat
!rjw            if(master) write(*,*) 'Warning: surf_map failed'
                write(*,*) 'Warning: surf_map_simple failed'
                q2(i,j) = 1.E8
                call mp_stop

            endif
         enddo
      enddo
      end subroutine map_to_cubed_simple

!#######################################################################
! reads the namelist file, write namelist to log file,
! and initializes constants

subroutine read_namelist

   integer :: unit, ierr, io
!   real    :: dtr, ght

!  read namelist

#ifdef INTERNAL_FILE_NML
   read  (input_nml_file, nml=surf_map_nml, iostat=io)
   ierr = check_nml_error(io,'surf_map_nml')
#else
   if ( file_exist('input.nml')) then
      unit = open_namelist_file ( )
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=surf_map_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'surf_map_nml')
      enddo
 10   call close_file (unit)
   endif
#endif

!  write version and namelist to log file

   if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=surf_map_nml)
   endif

end subroutine read_namelist


 end module fv_surf_map_mod
