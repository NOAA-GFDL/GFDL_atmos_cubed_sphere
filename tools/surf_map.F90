 module surf_map

      use constants_mod, only: grav
!     use tpcore,     only: copy_corners
      use grid_utils, only: great_circle_dist, latlon2xyz, v_prod,  &
                            sina_u, sina_v, g_sum, global_mx 
      use mp_mod,   only: domain, ng, is,js,ie,je, isd,jsd,ied,jed, &
                          mp_stop, mp_reduce_min, mp_reduce_max
      use mpp_domains_mod,   only : mpp_update_domains
      use fms_mod, only: file_exist, error_mesg, FATAL

      implicit none
      real pi
      private
      real, allocatable:: sgh_g(:,:), oro_g(:,:), zs_g(:,:)
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
! Position of edges of the box containing the original data point:

      integer          londim
      integer          latdim
      character*80:: oflnm
      character*80:: hgtflnm

      real dx1, dx2, dy1, dy2

      character*80     iflnm
      character*80     topoflnm
      real*4 fmin, fmax, vmax
      real*4, allocatable :: ftopo(:,:),  htopo(:,:)
      real, allocatable :: lon1(:),  lat1(:)
      integer i, j, n
      integer ncid, lonid, latid, ftopoid, htopoid
      integer status
      logical check_orig
      integer nlon, nlat
      real da_min, da_max, cd2, cd4, zmean
      integer fid

! Output the original 10 min NAVY data in grads readable format
      data             check_orig /.false./

      allocate ( oro_g(isd:ied, jsd:jed) )
      allocate ( sgh_g(isd:ied, jsd:jed) )
      allocate (  zs_g(is:ie, js:je) )

#ifdef MARS_GCM
!!!             Note that it is necessary to convert from km to m
!!!                     see ifdef MARS_GCM  below 
      fid = 22
      if(file_exist('INPUT/mars_topo')then
         open( Unit= fid, FILE= 'INPUT/mars_topo', FORM= 'unformatted')
         read( fid ) nlon, nlat

         if(master) write(*,*) 'Mars Terrain dataset dims=', nlon, nlat

         allocate ( htopo(nlon,nlat) )
         read( fid )  htopo

         if ( master ) then
              write(6,*) 'Check Hi-res Mars data ..'
!             write(*,*) 1.E3*htopo(1,1), 1.E3*htopo(nlon, nlat)
              fmax =  vmax(htopo,fmin,nlon,nlat,1)
              write(6,*) 'hmax=', fmax*1.E3
              write(6,*) 'hmin=', fmin*1.E3
         endif
         close(fid)
!     if(master) write(fid)  htopo


         allocate ( ftopo(nlon,nlat) )
         ftopo = 1.              ! all lands
      else
         call error_mesg ( 'surfdrv',  &
             'mars_topo not found in INPUT', FATAL )
      endif
#else
 
#ifndef USE_IEEE_DATA
      if ( npx > 65 ) then
           nlon = 10800
           nlat =  5400
           iflnm  = 'INPUT/topo2min.nc'
      else
           iflnm  = 'INPUT/topo5min.nc'
           nlon = 4320
           nlat = 2160
      endif

      if(master) write(*,*) 'USGS dataset = ', iflnm

      topoflnm = 'topo.bin'

! ... Open input netCDF file

      allocate ( ftopo(nlon,nlat) )
      allocate ( htopo(nlon,nlat) )

      if ( master ) write(*,*) 'Opening USGS datset file'

      status = nf_open (iflnm, NF_NOWRITE, ncid)
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
!         real*4 data (GrADS format)
!
      if (check_orig) then
        open (31, file=topoflnm, form='unformatted', &
            status='unknown', access='direct', recl=nlon*nlat*4)
        write (31, rec=1) ftopo
        write (31, rec=2) htopo
        close (31)
      endif
#else
      nlon = 10800
      nlat =  5400
      iflnm    = 'INPUT/topo2min.bin'

      if(master) write(*,*) 'USGS dataset (binary) = ', iflnm

      open (31, file=iflnm, form='unformatted', &
                status='unknown', access='direct', recl=nlon*nlat*4)

      allocate ( ftopo(nlon,nlat) )
      allocate ( htopo(nlon,nlat) )
      read (31, rec=1) ftopo
      read (31, rec=2) htopo
      close (31)
#endif
#endif

      allocate ( lat1(nlat+1) )
      allocate ( lon1(nlon+1) )

      pi = 4.0 * datan(1.0d0)

      dx1 = 2.*pi/real(nlon)
      dy1 = pi/real(nlat)

      do i=1,nlon+1
         lon1(i) = dx1 * (i-1)
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

      if ( npx>1000 ) then
! Use a quicker algorithm based on distance to cell center; this is a better approach
! if the input dataset resolution is close to the target resolution. Filtering is done
! by using a larger radius of inclusion.
! To be implemented:
      call map_to_cubed_dist(nlon, nlat, lat1, lon1, htopo, ftopo, grid, agrid,  &
                            phis, oro_g, sgh_g, master, npx, npy)
      else
      call map_to_cubed_raw(nlon, nlat, lat1, lon1, htopo, ftopo, grid, agrid,  &
                            phis, oro_g, sgh_g, master, npx, npy)
      endif

      if(master) write(*,*) 'map_to_cubed_raw: done'

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
         end do
      end do
!--------
! Filter:
!--------
      call global_mx(area, ng, da_min, da_max)

! diffusion coefficient:
      cd2 = 0.10 * da_min
      cd4 = 0.24 * da_min

      call global_mx(phis, ng, da_min, da_max)
      zmean = g_sum(zs_g, is, ie, js, je, ng, area, mode=1)

      if ( master )  &
           write(*,*) 'Before filter Phis min=', da_min, ' Max=', da_max,' Mean=',zmean

      if ( npx<=1000 ) then 
         call del4_cubed_sphere(npx, npy, phis, area, dx, dy, dxc, dyc, 2, cd4)
!        call del2_cubed_sphere(npx, npy, phis, area, dx, dy, dxc, dyc, 1, cd2)
      endif

      call global_mx(phis, ng, da_min, da_max)
      zmean = g_sum(phis(is:ie,js:je), is, ie, js, je, ng, area, mode=1)
      if ( master ) &
           write(*,*) 'After filter Phis min=', da_min, ' Max=', da_max, 'Mean=', zmean

      do j=js,je
         do i=is,ie
#ifdef MARS_GCM
            phis(i,j) = grav * 1.E3 * phis(i,j)  ! Convert km to meters
#else
            phis(i,j) =  grav * phis(i,j)
#endif
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

      if ( npx<=1000 ) then
         call del4_cubed_sphere(npx, npy, sgh_g, area, dx, dy, dxc, dyc, 2, cd4)
!        call del2_cubed_sphere(npx, npy, sgh_g, area, dx, dy, dxc, dyc, 1, cd2)
         call global_mx(sgh_g, ng, da_min, da_max)
         if ( master ) write(*,*) 'After filter SGH min=', da_min, ' Max=', da_max
         do j=js,je
            do i=is,ie
               sgh_g(i,j) = max(0., sgh_g(i,j))
            enddo
         enddo
      endif


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

      call mpp_update_domains(q,domain,whalo=1,ehalo=1,shalo=1,nhalo=1)

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

         if( n>1 ) call mpp_update_domains(q,domain,whalo=1,ehalo=1,shalo=1,nhalo=1)

         do j=js,je
            do i=is,ie+1
               ddx(i,j) = dy(i,j)*sina_u(i,j)*(q(i,j)-q(i-1,j))/dxc(i,j)
            enddo
         enddo

         do j=js,je+1
            do i=is,ie
               ddy(i,j) = dx(i,j)*sina_v(i,j)*(q(i,j)-q(i,j-1))/dyc(i,j)
            enddo
         enddo

         do j=js,je
            do i=is,ie
               q(i,j) = q(i,j) + cd/area(i,j)*(ddx(i+1,j)-ddx(i,j)+ddy(i,j+1)-ddy(i,j))
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

      call mpp_update_domains(q,domain,whalo=1,ehalo=1,shalo=1,nhalo=1)

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

  do n=1,nmax
     if( n/=1 ) call mpp_update_domains(q,domain,whalo=1,ehalo=1,shalo=1,nhalo=1)
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

      call mpp_update_domains(d2,domain,whalo=1,ehalo=1,shalo=1,nhalo=1)

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

      call mpp_update_domains(win,domain,whalo=1,ehalo=1,shalo=1,nhalo=1, complete=.false.)
      call mpp_update_domains(wou,domain,whalo=1,ehalo=1,shalo=1,nhalo=1, complete=.true.)

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




      subroutine map_to_cubed_dist(im, jm, lat1, lon1, q1, f1,  grid, agrid,  &
                                   q2, f2, h2, master, npx, npy)

! Input
      integer, intent(in):: im,jm         ! original dimensions
      integer, intent(in):: npx, npy
      logical, intent(in):: master
      real, intent(in):: lat1(jm+1)       ! original southern edge of the cell [-pi/2:pi/2]
      real, intent(in):: lon1(im+1)       ! original western edge of the cell [0:2*pi]
      real*4, intent(in):: q1(im,jm)      ! original data at center of the cell
      real*4, intent(in):: f1(im,jm)      !

      real, intent(in)::  grid(isd:ied+1, jsd:jed+1,2)
      real, intent(in):: agrid(isd:ied,   jsd:jed,  2)

! Output
      real, intent(out):: q2(isd:ied,jsd:jed) ! Mapped data at the target resolution
      real, intent(out):: f2(isd:ied,jsd:jed) ! oro
      real, intent(out):: h2(isd:ied,jsd:jed) ! variances of terrain

! Local
      real*4  qt(-im/32:im+im/32,jm)    ! ghosted east-west
      real*4  ft(-im/32:im+im/32,jm)    ! 
      real lon_g(-im/32:im+im/32)
      real lat_g(jm)

      real pc(3), p2(2), pp(3), grid3(3,is-ng:ie+ng+1, js-ng:je+ng+1)
      integer i,j, np
      integer ii, jj, i1, i2, j1, j2
      integer ifirst, ilast
      real ddeg, latitude, qsum, fsum, hsum, lon_w, lon_e, lat_s, lat_n, r2d
      real delg

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
            ft(i,j) = f1(i,j)
         enddo
         do i=-im/32,0
            qt(i,j) = qt(i+im,j)
            ft(i,j) = ft(i+im,j)
         enddo
         do i=im+1,im+im/32
            qt(i,j) = qt(i-im,j)
            ft(i,j) = ft(i-im,j)
         enddo
      enddo

      do j=js,je+1
         do i=is,ie+1
            call latlon2xyz(grid(i,j,1:2), grid3(1,i,j))
         enddo
      enddo

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
            fsum = 0.
            hsum = 0.
            call latlon2xyz(agrid(i,j,1:2), pc)

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
                  if (inside_p4(grid3(1,i,j), grid3(1,i+1,j), grid3(1,i+1,j+1), grid3(1,i,j+1), pc, pp)) then
                       np = np + 1
                       qsum = qsum + qt(ii,jj)
                       fsum = fsum + ft(ii,jj)
                       hsum = hsum + qt(ii,jj)**2
                  endif
               enddo
            enddo
! Compute weighted average:
            if ( np > 0 ) then
                 q2(i,j) = qsum / real(np)
                 f2(i,j) = fsum / real(np)
                 h2(i,j) = hsum / real(np) - q2(i,j)**2
            else                    ! the subdomain could be totally flat
                 if(master) write(*,*) 'Warning: surf_map failed'
                 q2(i,j) = 1.E8
                 call mp_stop
            endif
         enddo
      enddo  

      end subroutine map_to_cubed_dist




      subroutine map_to_cubed_raw(im, jm, lat1, lon1, q1, f1,  grid, agrid,  &
                                  q2, f2, h2, master, npx, npy)

! Input
      integer, intent(in):: im,jm         ! original dimensions
      integer, intent(in):: npx, npy
      logical, intent(in):: master
      real, intent(in):: lat1(jm+1)       ! original southern edge of the cell [-pi/2:pi/2]
      real, intent(in):: lon1(im+1)       ! original western edge of the cell [0:2*pi]
      real*4, intent(in):: q1(im,jm)      ! original data at center of the cell
      real*4, intent(in):: f1(im,jm)      !

      real, intent(in)::  grid(isd:ied+1, jsd:jed+1,2)
      real, intent(in):: agrid(isd:ied,   jsd:jed,  2)

! Output
      real, intent(out):: q2(isd:ied,jsd:jed) ! Mapped data at the target resolution
      real, intent(out):: f2(isd:ied,jsd:jed) ! oro
      real, intent(out):: h2(isd:ied,jsd:jed) ! variances of terrain

! Local
      real*4  qt(-im/32:im+im/32,jm)    ! ghosted east-west
      real*4  ft(-im/32:im+im/32,jm)    ! 
      real lon_g(-im/32:im+im/32)
      real lat_g(jm)

      real pc(3), p2(2), pp(3), grid3(3,is-ng:ie+ng+1, js-ng:je+ng+1)
      integer i,j, np
      integer ii, jj, i1, i2, j1, j2
      integer ifirst, ilast
      real ddeg, latitude, qsum, fsum, hsum, lon_w, lon_e, lat_s, lat_n, r2d
      real delg

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
            ft(i,j) = f1(i,j)
         enddo
         do i=-im/32,0
            qt(i,j) = qt(i+im,j)
            ft(i,j) = ft(i+im,j)
         enddo
         do i=im+1,im+im/32
            qt(i,j) = qt(i-im,j)
            ft(i,j) = ft(i-im,j)
         enddo
      enddo

      do j=js,je+1
         do i=is,ie+1
            call latlon2xyz(grid(i,j,1:2), grid3(1,i,j))
         enddo
      enddo

     if(master) write(*,*) 'surf_map: Search started ....'
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
            fsum = 0.
            hsum = 0.
            call latlon2xyz(agrid(i,j,1:2), pc)

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
                  if (inside_p4(grid3(1,i,j), grid3(1,i+1,j), grid3(1,i+1,j+1), grid3(1,i,j+1), pc, pp)) then
                       np = np + 1
                       qsum = qsum + qt(ii,jj)
                       fsum = fsum + ft(ii,jj)
                       hsum = hsum + qt(ii,jj)**2
                  endif
               enddo
            enddo
! Compute weighted average:
            if ( np > 0 ) then
                 q2(i,j) = qsum / real(np)
                 f2(i,j) = fsum / real(np)
                 h2(i,j) = hsum / real(np) - q2(i,j)**2
            else                    ! the subdomain could be totally flat
                 if(master) write(*,*) 'Warning: surf_map failed'
                 q2(i,j) = 1.E8
                 call mp_stop
            endif
         enddo
      enddo  

      end subroutine map_to_cubed_raw



      logical function inside_p4(p1, p2, p3, p4, pc, pp)
      real, intent(in):: p1(3), p2(3), p3(3), p4(3)
      real, intent(in):: pc(3), pp(3)
! Local:
      real v1(3), v2(3), vp(3)
      real a1, a2, aa
      integer k

! S-W:
      do k=1,3
         v1(k) = p2(k) - p1(k) 
         v2(k) = p4(k) - p1(k) 
         vp(k) = pp(k) - p1(k) 
      enddo

      aa = v_prod(v1, v2)
      a1 = v_prod(v1, vp)
      a2 = v_prod(v2, vp)

      if ( min(a1,a2) < aa ) then
           inside_p4 = .false.
           return
      endif

! N-E:
      do k=1,3
         v1(k) = p2(k) - p3(k) 
         v2(k) = p4(k) - p3(k) 
         vp(k) = pp(k) - p3(k) 
      enddo
      aa = v_prod(v1, v2)
      a1 = v_prod(v1, vp)
      a2 = v_prod(v2, vp)

      if ( min(a1,a2) < aa ) then
           inside_p4 = .false.
           return
      endif

!     inside_p4 = .true.
!     return
! Final extra check (to exclude points from opposite of the sphere)

      aa = (pp(1)-pc(1))**2 + (pp(2)-pc(2))**2 + (pp(3)-pc(3))**2
      a1 = (p1(1)-p3(1))**2 + (p1(2)-p3(2))**2 + (p1(3)-p3(3))**2
      a2 = (p2(1)-p4(1))**2 + (p2(2)-p4(2))**2 + (p2(3)-p4(3))**2

      if ( aa < max(a1, a2) ) then
           inside_p4 = .true.
      else
           inside_p4 = .false.
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
      real*4 pmin, pmax
      real*4 a(m,n,z)

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
      real*4, intent(in):: q1(im,jm)        ! original data at center of the cell
!rjw      real*4, intent(in):: f1(im,jm)        !

      real, intent(in)::  grid(is-ng:ie+ng+1, js-ng:je+ng+1,2)
      real, intent(in):: agrid(is-ng:ie+ng,   js-ng:je+ng,  2)

! Output
      real, intent(out):: q2(is-ng:ie+ng, js-ng:je+ng) ! Mapped data at the target resolution
!rjw      real, intent(out):: f2(isd:ied,jsd:jed) ! oro
!rjw      real, intent(out):: h2(isd:ied,jsd:jed) ! variances of terrain

! Local
      real*4  qt(-im/32:im+im/32,jm)    ! ghosted east-west
!rjw      real*4  ft(-im/32:im+im/32,jm)    ! 
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
            call latlon2xyz(agrid(i,j,1:2), pc)

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
                  if (inside_p4(grid3(1,i,j), grid3(1,i+1,j), grid3(1,i+1,j+1), grid3(1,i,j+1), pc, pp)) then
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

 end module surf_map
