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

!>@brief The module 'fv_tracker' contains the internal GFDL/NCEP vortex tracker
!adapted from HWRF internal vortex tracker, mainly based on the GFDL vortex
!tracker.

module fv_tracker_mod

#ifdef MOVING_NEST

   use constants_mod,       only: pi=>pi_8, rad_to_deg, deg_to_rad
   use time_manager_mod,    only: time_type, get_time, set_time, operator(+), &
                                  operator(-), operator(/), time_type_to_real
   use mpp_mod,             only: mpp_error, stdout, FATAL, WARNING, NOTE, &
                                  mpp_root_pe, mpp_npes, mpp_pe, mpp_chksum, &
                                  mpp_get_current_pelist, &
                                  mpp_set_current_pelist, mpp_sync
   use mpp_domains_mod,     only: mpp_get_data_domain, mpp_get_compute_domain
   use fv_arrays_mod,       only: fv_atmos_type, R_GRID
   use fv_diagnostics_mod,  only: fv_diag_init, fv_diag, fv_time, prt_maxmin, prt_height
   use fv_mp_mod,           only: is_master, &
                                  mp_reduce_sum, mp_reduce_max, mp_reduce_min, &
                                  mp_reduce_minval, mp_reduce_maxval, &
                                  mp_reduce_minloc, mp_reduce_maxloc

  implicit none
  private
  public :: fv_tracker_init, fv_tracker_center, fv_tracker_post_move, fv_tracker_move

  integer, parameter :: maxtp=11 ! number of tracker parameters

  real, parameter :: invE=0.36787944117 ! 1/e
  real, parameter :: searchrad_6=250.0 ! km - ignore data more than this far from domain center
  real, parameter :: searchrad_7=200.0 ! km - ignore data more than this far from domain center
  real, parameter :: uverrmax=225.0 ! For use in get_uv_guess
  real, parameter :: ecircum=40030.2 ! Earth's circumference (km) using erad=6371.e3
  real, parameter :: rads_vmag=120.0 ! max search radius for wind minimum
  real, parameter :: err_reg_init=300.0 ! max err at initial time (km)
  real, parameter :: err_reg_max=225.0 ! max err at other times (km)

  real, parameter :: errpmax=485.0 ! max stddev of track parameters
  real, parameter :: errpgro=1.25 ! stddev multiplier

  real, parameter :: max_wind_search_radius=searchrad_7 ! max radius for vmax search
  real, parameter :: min_mlsp_search_radius=searchrad_7 ! max radius for pmin search

  real, parameter :: km2nmi=0.539957, kn2mps=0.514444, mps2kn=1./kn2mps

contains

  subroutine fv_tracker_init(Atm)
    ! Initialize tracker variables in the Atm structure.
    implicit none
    type(fv_atmos_type), intent(inout) :: Atm

    call mpp_error(NOTE, 'fv_tracker_init')

    Atm%track_stderr_m1=-99.9
    Atm%track_stderr_m2=-99.9
    Atm%track_stderr_m3=-99.9
  ! Atm%track_n_old=0
  ! Atm%track_old_lon=0
  ! Atm%track_old_lat=0
  ! Atm%track_old_ntsd=0

    Atm%tracker_angle=0
    Atm%tracker_fixlon=-999.0
    Atm%tracker_fixlat=-999.0
    Atm%tracker_ifix=-99
    Atm%tracker_jfix=-99
    Atm%tracker_havefix=.false.
    Atm%tracker_gave_up=.false.
    Atm%tracker_pmin=-99999.
    Atm%tracker_vmax=-99.
    Atm%tracker_rmw=-99.

    Atm%track_have_guess=.false.
    Atm%track_guess_lat=-999.0
    Atm%track_guess_lon=-999.0

  end subroutine fv_tracker_init

  subroutine fv_tracker_center(Atm)
    ! Top-level entry to the internal GFDL/NCEP vortex tracker. Finds the center of
    ! the storm in the specified Atm and updates the Atm variables.
    ! Will do nothing and return immediately if
    ! Atm%tracker_gave_up=.true.
    implicit none
    type(fv_atmos_type), intent(inout) :: Atm

    integer :: ids,ide,jds,jde,kds,kde
    integer :: ims,ime,jms,jme,kms,kme
    integer :: ips,ipe,jps,jpe,kps,kpe

    call mpp_error(NOTE, 'fv_tracker_center')
    CALL get_ijk_from_domain(Atm,         &
         ids, ide, jds, jde, kds, kde,    &
         ims, ime, jms, jme, kms, kme,    &
         ips, ipe, jps, jpe, kps, kpe    )

    call ntc_impl(Atm,                    &
         ids, ide, jds, jde, kds, kde,    &
         ims, ime, jms, jme, kms, kme,    &
         ips, ipe, jps, jpe, kps, kpe    )

  end subroutine fv_tracker_center

  subroutine ntc_impl(Atm,      &
       ids,ide,jds,jde,kds,kde, &
       ims,ime,jms,jme,kms,kme, &
       ips,ipe,jps,jpe,kps,kpe)
    ! This is the main entry point to the tracker.  It is most similar
    ! to the function "tracker" in the GFDL/NCEP vortex tracker.

    implicit none
    type(fv_atmos_type), intent(inout) :: Atm
    integer, intent(in) :: ids,ide,jds,jde,kds,kde
    integer, intent(in) :: ims,ime,jms,jme,kms,kme
    integer, intent(in) :: ips,ipe,jps,jpe,kps,kpe

    real :: dxdymean, sumdxa, sumdya
    integer :: i, j, iweights, ip

    integer :: iguess, jguess ! first guess location
    real :: latguess, longuess ! same, but in lat & lon

    integer :: iuvguess, juvguess ! "second guess" location using everything except wind maxima
    real :: srsq
    integer :: ifinal, jfinal
    real :: latfinal, lonfinal
    integer :: ierr
    integer :: icen(maxtp), jcen(maxtp) ! center locations for each parameter
    real :: loncen(maxtp), latcen(maxtp) ! lat, lon locations in degrees
    logical :: calcparm(maxtp) ! do we have a valid center location for this parameter?
    real :: max_wind, min_pres ! for ATCF output
    real :: rcen(maxtp) ! center value (max wind, min mslp, etc.)
    character*255 :: message
    logical :: north_hemi ! true = northern hemisphere
    logical :: have_guess ! first guess is available
    real :: guessdist, guessdeg ! first guess distance to nearest point on grid
    real :: latnear, lonnear ! nearest point in grid to first guess

    ! icen,jcen: Same meaning as clon, clat in tracker, but uses i and
    ! j indexes of the center instead of lat/lon.  Tracker comment:
    !            Holds the coordinates for the center positions for
    !            all storms at all times for all parameters.
    !            (max_#_storms, max_fcst_times, max_#_parms).
    !            For the third position (max_#_parms), here they are:
    !             1: Relative vorticity at 850 mb
    !             2: Relative vorticity at 700 mb
    !             3: Vector wind magnitude at 850 mb
    !             4: NOT CURRENTLY USED
    !             5: Vector wind magnitude at 700 mb
    !             6: NOT CURRENTLY USED
    !             7: Geopotential height at 850 mb
    !             8: Geopotential height at 700 mb
    !             9: Mean Sea Level Pressure
    !            10: Vector wind magnitude at 10 m
    !            11: Relative vorticity at 10 m

    call mpp_error(NOTE, 'ntc_impl')

    ! Initialize center information to invalid values for all centers:
    icen=-99
    jcen=-99
    latcen=9e9
    loncen=9e9
    rcen=9e9
    calcparm=.false.
    if(Atm%neststruct%vortex_tracker==6) then
       srsq=searchrad_6*searchrad_6*1e6
    else
       srsq=searchrad_7*searchrad_7*1e6
    endif

    ! Estimate the domain wide mean grid spacing in km
    sumdxa=0.0
    sumdya=0.0
    do j=jps,min(jde-1,jpe)
       do i=ips,min(ide-1,ipe)
          sumdxa=sumdxa+Atm%gridstruct%dxa(i,j)
          sumdya=sumdya+Atm%gridstruct%dya(i,j)
       enddo
    enddo
    call mp_reduce_sum(sumdxa)
    call mp_reduce_sum(sumdya)
    dxdymean=0.5*(sumdxa + sumdya)/((ide-ids) * (jde-jds)) / 1000.0

    ! Get the square of the approximate distance to the domain center
    ! at all points:
    call get_distsq(Atm, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe)

    ! Get the first guess from the prior nest motion timestep:
    have_guess=Atm%track_have_guess
    if(have_guess) then
       ! We have a first guess center.  We have to translate it to gridpoint space.
       longuess=Atm%track_guess_lon
       latguess=Atm%track_guess_lat
       call get_nearest_lonlat(Atm,iguess,jguess,ierr,longuess,latguess, &
            ids,ide, jds,jde, kds,kde, &
            ims,ime, jms,jme, kms,kme, &
            ips,ipe, jps,jpe, kps,kpe, &
            lonnear, latnear)
       if(ierr==0) then
          call calcdist(longuess,latguess, lonnear,latnear, guessdist,guessdeg)
          if(guessdist > Atm%neststruct%refinement*dxdymean) then
108          format('WARNING: guess lon=',F0.3,',lat=',F0.3, &
                  ' too far (',F0.3,'km) from nearest point lon=',F0.3,',lat=',F0.3, &
                  '.  Will use domain center as first guess.')
             write(message,108) Atm%track_guess_lon,Atm%track_guess_lat, &
                  guessdist,lonnear,latnear
             call mpp_error(NOTE, message)
             have_guess=.false. ! indicate that the first guess is unusable
          else
             latguess=latnear
             longuess=lonnear
          endif
       else
          have_guess=.false. ! indicate that the first guess is unusable.
109       format('WARNING: guess lon=',F0.3,',lat=',F0.3, &
                  ' does not exist in this domain.  Will use domain center as first guess.')
          write(message,109) Atm%track_guess_lon,Atm%track_guess_lat
          call mpp_error(NOTE, message)
       endif
   endif

   ! If we could not get the first guess from the prior nest motion
   ! timestep, then use the default first guess: the domain center.
    if(Atm%neststruct%vortex_tracker==6 .or. .not.have_guess) then
       ! vt=6: hard coded first-guess center is domain center:
       ! vt=7: first guess comes from prior timestep
       !     Initial first guess is domain center.
       !     Backup first guess is domain center if first guess is unusable.
       iguess=(ide-ids)/2+ids
       jguess=(jde-jds)/2+jds
       if(Atm%neststruct%vortex_tracker==7) then
          call mpp_error(NOTE, 'Using domain center as first guess since no valid first guess is available.')
       endif
       call get_lonlat(Atm,iguess,jguess,longuess,latguess,ierr, &
            ids,ide, jds,jde, kds,kde, &
            ims,ime, jms,jme, kms,kme, &
            ips,ipe, jps,jpe, kps,kpe)
       if(ierr/=0) then
          call mpp_error(FATAL, "ERROR: center of domain is not inside the domain")
       endif
       have_guess=.true.
    endif

    if(.not.have_guess) then
       call mpp_error(FATAL, "INTERNAL ERROR: No first guess is available (should never happen).")
    endif

    north_hemi = latguess>0.0

    ! Find the centers of all fields except the wind minima:
    call find_center(Atm,Atm%vort850,srsq, &
         icen(1),jcen(1),rcen(1),calcparm(1),loncen(1),latcen(1),dxdymean,'zeta', &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe, north_hemi=north_hemi)
    call find_center(Atm,Atm%vort700,srsq, &
         icen(2),jcen(2),rcen(2),calcparm(2),loncen(2),latcen(2),dxdymean,'zeta', &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe, north_hemi=north_hemi)
    call find_center(Atm,Atm%z850,srsq, &
         icen(7),jcen(7),rcen(7),calcparm(7),loncen(7),latcen(7),dxdymean,'hgt', &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe)
    call find_center(Atm,Atm%z700,srsq, &
         icen(8),jcen(8),rcen(8),calcparm(8),loncen(8),latcen(8),dxdymean,'hgt', &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe)
    call find_center(Atm,Atm%slp,srsq, &
         icen(9),jcen(9),rcen(9),calcparm(9),loncen(9),latcen(9),dxdymean,'slp', &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe)
    call find_center(Atm,Atm%vort10m,srsq, &
         icen(11),jcen(11),rcen(11),calcparm(11),loncen(11),latcen(11),dxdymean,'zeta', &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe, north_hemi=north_hemi)

    ! Get a guess center location for the wind minimum searches:
    call get_uv_guess(Atm,icen,jcen,loncen,latcen,calcparm, &
         iguess,jguess,longuess,latguess,iuvguess,juvguess, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe)

    ! Find wind minima.  Requires a first guess center:
    windmin: if(Atm%neststruct%vortex_tracker==6) then
       call find_center(Atm,Atm%spd850,srsq, &
            icen(3),jcen(3),rcen(3),calcparm(3),loncen(3),latcen(3),dxdymean,'wind', &
            ids,ide,jds,jde,kds,kde, &
            ims,ime,jms,jme,kms,kme, &
            ips,ipe,jps,jpe,kps,kpe, &
            iuvguess=iuvguess, juvguess=juvguess)
       call find_center(Atm,Atm%spd700,srsq, &
            icen(5),jcen(5),rcen(5),calcparm(5),loncen(5),latcen(5),dxdymean,'wind', &
            ids,ide,jds,jde,kds,kde, &
            ims,ime,jms,jme,kms,kme, &
            ips,ipe,jps,jpe,kps,kpe, &
            iuvguess=iuvguess, juvguess=juvguess)
       call find_center(Atm,Atm%spd10m,srsq, &
            icen(10),jcen(10),rcen(10),calcparm(10),loncen(10),latcen(10),dxdymean,'wind', &
            ids,ide,jds,jde,kds,kde, &
            ims,ime,jms,jme,kms,kme, &
            ips,ipe,jps,jpe,kps,kpe, &
            iuvguess=iuvguess, juvguess=juvguess)
    else
       call get_uv_center(Atm,Atm%spd850, &
            icen(3),jcen(3),rcen(3),calcparm(3),loncen(3),latcen(3),dxdymean,'wind', &
            ids,ide,jds,jde,kds,kde, &
            ims,ime,jms,jme,kms,kme, &
            ips,ipe,jps,jpe,kps,kpe, &
            iuvguess=iuvguess, juvguess=juvguess)
       call get_uv_center(Atm,Atm%spd700, &
            icen(5),jcen(5),rcen(5),calcparm(5),loncen(5),latcen(5),dxdymean,'wind', &
            ids,ide,jds,jde,kds,kde, &
            ims,ime,jms,jme,kms,kme, &
            ips,ipe,jps,jpe,kps,kpe, &
            iuvguess=iuvguess, juvguess=juvguess)
       call get_uv_center(Atm,Atm%spd10m, &
            icen(10),jcen(10),rcen(10),calcparm(10),loncen(10),latcen(10),dxdymean,'wind', &
            ids,ide,jds,jde,kds,kde, &
            ims,ime,jms,jme,kms,kme, &
            ips,ipe,jps,jpe,kps,kpe, &
            iuvguess=iuvguess, juvguess=juvguess)
    endif windmin

    ! Get a final guess center location:
    call fixcenter(Atm,icen,jcen,calcparm,loncen,latcen, &
         iguess,jguess,longuess,latguess, &
         ifinal,jfinal,lonfinal,latfinal, &
         north_hemi, &
         ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         ips,ipe, jps,jpe, kps,kpe)

    Atm%tracker_fixes=0
    do ip=1,maxtp
       if(calcparm(ip)) then
300       format('Parameter ',I0,': i=',I0,' j=',I0,' lon=',F0.2,' lat=',F0.2)
          !write(0,300) ip,icen(ip),jcen(ip),loncen(ip),latcen(ip)
          if(icen(ip)>=ips .and. icen(ip)<=ipe &
               .and. jcen(ip)>=jps .and. jcen(ip)<=jpe) then
             Atm%tracker_fixes(icen(ip),jcen(ip))=ip
          endif
       else
301       format('Parameter ',I0,' invalid')
          !write(0,301) ip
       endif
    enddo

    if(iguess>=ips .and. iguess<=ipe .and. jguess>=jps .and. jguess<=jpe) then
       Atm%tracker_fixes(iguess,jguess)=-1
201    format('First guess: i=',I0,' j=',I0,' lon=',F0.2,' lat=',F0.2)
       !write(0,201) iguess,jguess,longuess,latguess
    endif

    if(iuvguess>=ips .and. iuvguess<=ipe .and. juvguess>=jps .and. juvguess<=jpe) then
       Atm%tracker_fixes(iuvguess,juvguess)=-2
202    format('UV guess: i=',I0,' j=',I0)
       !write(0,202) iguess,jguess
    endif

1000 format('Back with final lat/lon at i=',I0,' j=',I0,' lon=',F0.3,' lat=',F0.3)
    !write(0,1000) ifinal,jfinal,lonfinal,latfinal

    if(ifinal>=ips .and. ifinal<=ipe .and. jfinal>=jps .and. jfinal<=jpe) then
       Atm%tracker_fixes(ifinal,jfinal)=-3
203    format('Final fix: i=',I0,' j=',I0,' lon=',F0.2,' lat=',F0.2)
       !write(0,201) ifinal,jfinal,lonfinal,latfinal
    endif

    call get_tracker_distsq(Atm, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe)

    call get_wind_pres_intensity(Atm, &
         Atm%tracker_pmin,Atm%tracker_vmax,Atm%tracker_rmw, &
         max_wind_search_radius, min_mlsp_search_radius, &
         lonfinal,latfinal, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe)

    if(is_master()) then
       call output_partial_atcfunix(Atm, &
            ids,ide,jds,jde,kds,kde, &
            ims,ime,jms,jme,kms,kme, &
            ips,ipe,jps,jpe,kps,kpe)
    endif
  end subroutine ntc_impl

  subroutine get_ijk_from_domain(Atm,  &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         ips, ipe, jps, jpe, kps, kpe )
    ! This is
    implicit none
    type(fv_atmos_type), intent(in) :: Atm

    integer, intent(out) :: ids,ide,jds,jde,kds,kde
    integer, intent(out) :: ims,ime,jms,jme,kms,kme
    integer, intent(out) :: ips,ipe,jps,jpe,kps,kpe

    ids = 1
    ide = Atm%npx
    jds = 1
    jde = Atm%npy
    kds = 1
    kde = Atm%npz
    call mpp_get_data_domain(Atm%domain, ims, ime, jms, jme)
    kms = 1
    kme = Atm%npz
    call mpp_get_compute_domain(Atm%domain, ips, ipe, jps, jpe)
    kps = 1
    kpe = Atm%npz
  end subroutine get_ijk_from_domain

  subroutine get_nearest_lonlat(Atm,iloc,jloc,ierr,lon,lat, &
               ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme, &
               ips,ipe, jps,jpe, kps,kpe, &
               lonnear, latnear)
    ! Finds the nearest point in the domain to the specified lon,lat
    ! location.
    implicit none
    type(fv_atmos_type), intent(inout) :: Atm
    integer, intent(in) :: ids,ide,jds,jde,kds,kde
    integer, intent(in) :: ims,ime,jms,jme,kms,kme
    integer, intent(in) :: ips,ipe,jps,jpe,kps,kpe
    integer, intent(out) :: iloc,jloc,ierr
    real, intent(in) :: lon,lat
    real :: dx,dy,d,dmin, zdummy, latmin,lonmin
    integer :: i,j,imin,jmin
    real, intent(out), optional :: latnear, lonnear

    zdummy=42
    dmin=9e9
    imin=-99
    jmin=-99
    latmin=9e9
    lonmin=9e9
    ierr=0
    do j=jps,min(jde-1,jpe)
       do i=ips,min(ide-1,ipe)
          dy=abs(lat-Atm%gridstruct%agrid(i,j,2)*rad_to_deg)
          dx=abs(mod(3600.+180.+(lon-Atm%gridstruct%agrid(i,j,1)*rad_to_deg),360.)-180.)
          d=dx*dx+dy*dy
          if(d<dmin) then
             dmin=d
             imin=i
             jmin=j
             latmin=Atm%gridstruct%agrid(i,j,2)*rad_to_deg
             lonmin=Atm%gridstruct%agrid(i,j,1)*rad_to_deg
          endif
       enddo
    enddo

    call mp_reduce_minloc(dmin,latmin,lonmin,zdummy,imin,jmin)
    if(imin<0 .or. jmin<0) then
       ierr=-99
    else
       iloc=imin ; jloc=jmin
    endif
    if(present(latnear)) latnear=latmin
    if(present(lonnear)) lonnear=lonmin
  end subroutine get_nearest_lonlat

  subroutine output_partial_atcfunix(Atm, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         its,ite,jts,jte,kts,kte)
    ! This outputs to a format that can be easily converted to ATCF,
    ! using units used by ATCF.
    implicit none
    type(fv_atmos_type), intent(inout) :: Atm
    integer, intent(in) :: ids,ide,jds,jde,kds,kde
    integer, intent(in) :: ims,ime,jms,jme,kms,kme
    integer, intent(in) :: its,ite,jts,jte,kts,kte
    real :: sec
    character*255 message

    sec=time_type_to_real(Atm%Time-Atm%Time_Init)
313 format(F11.2,", ",                                  &
           "W10 = ",F7.3," kn, PMIN = ",F8.3," mbar, ", &
           "LAT = ",F6.3,A1,", LON = ",F7.3,A1,", ",    &
           "RMW = ",F7.3," nmi")
    write(Atm%neststruct%outatcf_lun+Atm%grid_number,313) sec,   &
         Atm%tracker_vmax*mps2kn,Atm%tracker_pmin/100.,          &
         abs(Atm%tracker_fixlat),get_lat_ns(Atm%tracker_fixlat), &
         abs(Atm%tracker_fixlon),get_lon_ew(Atm%tracker_fixlon), &
         Atm%tracker_rmw*km2nmi
    ! write(message,313) sec,                                      &
    !      Atm%tracker_vmax*mps2kn,Atm%tracker_pmin/100.,          &
    !      abs(Atm%tracker_fixlat),get_lat_ns(Atm%tracker_fixlat), &
    !      abs(Atm%tracker_fixlon),get_lon_ew(Atm%tracker_fixlon), &
    !      Atm%tracker_rmw*km2nmi
    ! call mpp_error(NOTE, message)
  end subroutine output_partial_atcfunix

  subroutine get_wind_pres_intensity(Atm, &
       min_mslp,max_wind,rmw, &
       max_wind_search_radius, min_mlsp_search_radius, clon,clat, &
       ids,ide,jds,jde,kds,kde, &
       ims,ime,jms,jme,kms,kme, &
       its,ite,jts,jte,kts,kte)
    ! This determines the maximum wind, RMW and minimum mslp in the domain.
    implicit none
    type(fv_atmos_type), intent(inout) :: Atm
    real, intent(out) :: min_mslp,max_wind,rmw
    real, intent(in) :: max_wind_search_radius, min_mlsp_search_radius,clon,clat
    integer, intent(in) :: ids,ide,jds,jde,kds,kde
    integer, intent(in) :: ims,ime,jms,jme,kms,kme
    integer, intent(in) :: its,ite,jts,jte,kts,kte

    real :: localextreme,globalextreme, sdistsq,windsq
    real :: globallat,globallon,degrees
    integer :: locali,localj,globali,globalj,ierr,i,j

    ! Get the MSLP minimum location and determine if what we found is
    ! still a storm:
    localextreme=9e9
    locali=-1
    localj=-1
    sdistsq=min_mlsp_search_radius*min_mlsp_search_radius*1e6
    do j=jts,min(jte,jde-1)
       do i=its,min(ite,ide-1)
          if(Atm%slp(i,j)<localextreme .and. &
               Atm%tracker_distsq(i,j)<sdistsq) then
             localextreme=Atm%slp(i,j)
             locali=i
             localj=j
          endif
       enddo
    enddo

    globalextreme=localextreme
    globali=locali
    globalj=localj
    call mp_reduce_minval(globalextreme,globali,globalj)
    min_mslp=globalextreme
    if(globali<0 .or. globalj<0) then
       call mpp_error(WARNING, "WARNING: No mslp values found that were less than 9*10^9.")
       min_mslp=-999
    endif

    ! Get the wind maximum location.  Note that we're using the square
    ! of the wind until after the loop to avoid the sqrt() call.
    localextreme=-9e9
    locali=-1
    localj=-1
    sdistsq=max_wind_search_radius*max_wind_search_radius*1e6
    do j=jts,min(jte,jde-1)
       do i=its,min(ite,ide-1)
          if(Atm%tracker_distsq(i,j)<sdistsq) then
             windsq=Atm%u10m(i,j)*Atm%u10m(i,j) + &
                    Atm%v10m(i,j)*Atm%v10m(i,j)
             if(windsq>localextreme) then
                localextreme=windsq
                locali=i
                localj=j
             endif
          endif
       enddo
    enddo
    if(localextreme>0) localextreme=sqrt(localextreme)

    globalextreme=localextreme
    globali=locali
    globalj=localj
    call mp_reduce_maxval(globalextreme,globali,globalj)

    call get_lonlat(Atm,globali,globalj,globallon,globallat,ierr, &
         ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         its,ite, jts,jte, kts,kte)
    if(ierr/=0) then
       call mpp_error(WARNING, "WARNING: Unable to find location of wind maximum.")
       rmw=-99
    else
       call calcdist(clon,clat,globallon,globallat,rmw,degrees)
    end if

    ! Get the guess location for the next time:
    max_wind=globalextreme
    if(globali<0 .or. globalj<0) then
       call mpp_error(WARNING, "WARNING: No wind values found that were greater than -9*10^9.")
       min_mslp=-999
    endif

  end subroutine get_wind_pres_intensity

  subroutine fixcenter(Atm,icen,jcen,calcparm,loncen,latcen, &
       iguess,jguess,longuess,latguess, &
       ifinal,jfinal,lonfinal,latfinal, &
       north_hemi, &
       ids,ide, jds,jde, kds,kde, &
       ims,ime, jms,jme, kms,kme, &
       ips,ipe, jps,jpe, kps,kpe)
    ! This is the same as "fixcenter" in gettrk_main.  Original comment:
    !
    ! ABSTRACT: This subroutine loops through the different parameters
    !           for the input storm number (ist) and calculates the
    !           center position of the storm by taking an average of
    !           the center positions obtained for those parameters.
    !           First we check to see which parameters are within a
    !           max error range (errmax), and we discard those that are
    !           not within that range.  Of the remaining parms, we get
    !           a mean position, and then we re-calculate the position
    !           by giving more weight to those estimates that are closer
    !           to this mean first-guess position estimate.

    ! Arguments: Input:
    ! grid - the grid being processed
    ! icen,jcen - arrays of center gridpoint locations
    ! calcperm - array of center validity flags (true = center is valid)
    ! loncen,latcen - center geographic locations
    ! iguess,jguess - first guess gridpoint location
    ! longuess,latguess - first guess geographic location

    ! Arguments: Output:
    ! ifinal,jfinal - final center gridpoint location
    ! lonfinal,latfinal - final center geographic location

    ! Arguments: Optional input:
    ! north_hemi - true = northern hemisphere, false=south

    implicit none
    integer, intent(in) :: &
         ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         ips,ipe, jps,jpe, kps,kpe
    type(fv_atmos_type), intent(inout) :: Atm
    integer, intent(in) :: icen(maxtp), jcen(maxtp)
    real, intent(in) :: loncen(maxtp), latcen(maxtp)
    logical, intent(inout) :: calcparm(maxtp)

    integer, intent(in) :: iguess,jguess
    real, intent(in) :: latguess,longuess

    integer, intent(inout) :: ifinal,jfinal
    real, intent(inout) :: lonfinal,latfinal

    logical, intent(in), optional :: north_hemi

    character*255 :: message
    real :: errdist(maxtp),avgerr,errmax,errinit,xavg_stderr
    real :: dist,degrees, total
    real :: minutes,hours,trkerr_avg,dist_from_mean(maxtp),wsum
    integer :: ip,itot4next,iclose,count,ifound,ierr
    integer(kind=8) :: isum,jsum
    real :: irsum,jrsum,errtmp,devia,wtpos
    real :: xmn_dist_from_mean, stderr_close
    logical use4next(maxtp)

    ! Determine forecast hour:
    hours=time_type_to_real(Atm%Time-Atm%Time_Init)/3600.

    ! Decide maximum values for distance and std. dev.:
    if(hours<0.5) then
       errmax=err_reg_init
       errinit=err_reg_init
    else
       errmax=err_reg_max
       errinit=err_reg_max
    endif

    if(hours>4.) then
       xavg_stderr = ( Atm%track_stderr_m1 + &
            Atm%track_stderr_m2 + Atm%track_stderr_m3 ) / 3.0
    elseif(hours>3.) then
       xavg_stderr = ( Atm%track_stderr_m1 + Atm%track_stderr_m2 ) / 2.0
    elseif(hours>2.) then
       xavg_stderr = Atm%track_stderr_m1
    endif

    if(hours>2.) then
       errtmp = 3.0*xavg_stderr*errpgro
       errmax = max(errtmp,errinit)
       errtmp = errpmax
       errmax = min(errmax,errtmp)
    endif

    ! Initialize loop variables:
    errdist=0.0
    use4next=.false.
    trkerr_avg=0
    itot4next=0
    iclose=0
    isum=0
    jsum=0
    ifound=0

    !write(0,*) 'errpmax=',errpmax
    !write(0,*) 'errmax=',errmax

500 format('Parm ip=',I0,' dist=',F0.3)
501 format('  too far, but discard')
    do ip=1,maxtp
       if(ip==4 .or. ip==6) then
          calcparm(ip)=.false.
          cycle
       elseif(calcparm(ip)) then
          ifound=ifound+1
          call calcdist(longuess,latguess,loncen(ip),latcen(ip),dist,degrees)
          errdist(ip)=dist
          !write(0,500) ip,dist
          if(dist<=errpmax) then
             if(ip==3 .or. ip==5 .or. ip==10) then
                use4next(ip)=.false.
                !write(0,'(A)') '  within range but discard: errpmax'
             else
                !write(0,'(A)') '  within range and keep: errpmax'
                use4next(ip)=.true.
                trkerr_avg=trkerr_avg+dist
                itot4next=itot4next+1
             endif
          endif
          if(dist<=errmax) then
502          format('  apply i=',I0,' j=',I0)
             !write(0,502) icen(ip),jcen(ip)
             iclose=iclose+1
             isum=isum+icen(ip)
             jsum=jsum+jcen(ip)
503          format(' added things isum=',I0,' jsum=',I0,' iclose=',I0)
             !write(0,503) isum,jsum,iclose
          else
             !write(0,*) '  discard; too far: errmax'
             calcparm(ip)=.false.
          endif
       endif
    enddo

    if(ifound<=0) then
       call mpp_error(NOTE, 'The tracker could not find the centers for any parameters. &
                             Thus, a center position could not be obtained for this storm.')
       goto 999
    endif

    if(iclose<=0) then
200    format('No storms are within errmax=',F0.1,'km of the parameters')
       write(message,200) errmax
       call mpp_error(NOTE, message)
       goto 999
    endif

    ifinal=real(isum)/real(iclose)
    jfinal=real(jsum)/real(iclose)

504 format(' calculated ifinal, jfinal: ifinal=',I0,' jfinal=',I0,' isum=',I0,' jsum=',I0,' iclose=',I0)
    !write(0,504) ifinal,jfinal,isum,jsum,iclose

    call get_lonlat(Atm,ifinal,jfinal,lonfinal,latfinal,ierr, &
         ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         ips,ipe, jps,jpe, kps,kpe)
    if(ierr/=0) then
       call mpp_error(NOTE, 'Gave up on finding the storm location due to error in get_lonlat (1).')
       goto 999
    endif

    count=0
    dist_from_mean=0.0
    total=0.0
    do ip=1,maxtp
       if(calcparm(ip)) then
          call calcdist(lonfinal,latfinal,loncen(ip),latcen(ip),dist,degrees)
          dist_from_mean(ip)=dist
          total=total+dist
          count=count+1
       endif
    enddo
    xmn_dist_from_mean=total/real(count)

    do ip=1,maxtp
       if(calcparm(ip)) then
          total=total+(xmn_dist_from_mean-dist_from_mean(ip))**2
       endif
    enddo
    if(count<2) then
       stderr_close=0.0
    else
       stderr_close=max(1.0,sqrt(1./(count-1) * total))
    endif

    if(calcparm(1) .or. calcparm(2) .or. calcparm(7) .or. &
         calcparm(8) .or. calcparm(9) .or. calcparm(11)) then
       continue
    else
       ! Message copied straight from tracker:
       call mpp_error(NOTE, 'In fixcenter, STOPPING PROCESSING for this storm.  The reason is that')
       call mpp_error(NOTE, 'none of the fix locations for parms z850, z700, zeta 850, zeta 700')
       call mpp_error(NOTE, 'MSLP or sfc zeta were within a reasonable distance of the guess location.')
       goto 999
    endif

    ! Recalculate the final center location using weights
    if(stderr_close<5.0) then
       ! Old code forced a minimum of 5.0 stddev
       stderr_close=5.0
    endif
    irsum=0
    jrsum=0
    wsum=0
    do ip=1,maxtp
       if(calcparm(ip)) then
          devia=max(1.0,dist_from_mean(ip)/stderr_close)
          wtpos=exp(-devia/3.)
          irsum=icen(ip)*wtpos+irsum
          jrsum=jcen(ip)*wtpos+jrsum
          wsum=wtpos+wsum
1100      format(' Adding parm: devia=',F0.3,' wtpos=',F0.3,' irsum=',F0.3,' jrsum=',F0.3,' wsum=',F0.3)
          !write(0,1100) devia,wtpos,irsum,jrsum,wsum
       endif
    enddo
    ifinal=nint(real(irsum)/real(wsum))
    jfinal=nint(real(jrsum)/real(wsum))
    call get_lonlat(Atm,ifinal,jfinal,lonfinal,latfinal,ierr, &
         ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         ips,ipe, jps,jpe, kps,kpe)
    if(ierr/=0) then
       call mpp_error(NOTE, 'Gave up on finding the storm location due to error in get_lonlat (2).')
       goto 999
    endif

    ! Store the lat/lon location:
    Atm%tracker_fixlon=lonfinal
    Atm%tracker_fixlat=latfinal
    Atm%tracker_ifix=ifinal
    Atm%tracker_jfix=jfinal
    Atm%tracker_havefix=.true.

1000 format('Stored lat/lon at i=',I0,' j=',I0,' lon=',F0.3,' lat=',F0.3)
    !write(0,1000) ifinal,jfinal,lonfinal,latfinal


    if(nint(hours) > Atm%track_last_hour ) then
       ! It is time to recalculate the std. dev. of the track:
       count=0
       dist_from_mean=0.0
       total=0.0
       do ip=1,maxtp
          if(calcparm(ip)) then
             call calcdist(lonfinal,latfinal,loncen(ip),loncen(ip),dist,degrees)
             dist_from_mean(ip)=dist
             total=total+dist
             count=count+1
          endif
       enddo
       xmn_dist_from_mean=total/real(count)

       do ip=1,maxtp
          if(calcparm(ip)) then
             total=total+(xmn_dist_from_mean-dist_from_mean(ip))**2
          endif
       enddo
       if(count<2) then
          stderr_close=0.0
       else
          stderr_close=max(1.0,sqrt(1./(count-1) * total))
       endif

       Atm%track_stderr_m3=Atm%track_stderr_m2
       Atm%track_stderr_m2=Atm%track_stderr_m1
       Atm%track_stderr_m1=stderr_close
       Atm%track_last_hour=nint(hours)
    endif

    !write(0,*) 'got to return'
    return

    ! We jump here if we're giving up on finding the center
999 continue
    ! Use domain center as storm location
    Atm%tracker_ifix=(ide-ids)/2+ids
    Atm%tracker_jfix=(jde-jds)/2+jds
    Atm%tracker_havefix=.false.
    Atm%tracker_gave_up=.true.
    call get_lonlat(Atm,ifinal,jfinal,lonfinal,latfinal,ierr, &
         ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         ips,ipe, jps,jpe, kps,kpe)
    if(ierr/=0) then
       call mpp_error(FATAL, 'Center of domain is not in domain (!?)')
       goto 999
    endif

    Atm%tracker_fixlon=-999.0
    Atm%tracker_fixlat=-999.0

  end subroutine fixcenter

  subroutine get_uv_guess(Atm,icen,jcen,loncen,latcen,calcparm, &
       iguess,jguess,longuess,latguess,iout,jout, &
       ids,ide,jds,jde,kds,kde, &
       ims,ime,jms,jme,kms,kme, &
       its,ite,jts,jte,kts,kte)
    ! This is a rewrite of the gettrk_main.f get_uv_guess.  Original comment:
    ! ABSTRACT: The purpose of this subroutine is to get a modified
    !           first guess lat/lon position before searching for the
    !           minimum in the wind field.  The reason for doing this is
    !           to better refine the guess and avoid picking up a wind
    !           wind minimum far away from the center.  So, use the
    !           first guess position (and give it strong weighting), and
    !           then also use the  fix positions for the current time
    !           (give the vorticity centers stronger weighting as well),
    !           and then take the average of these positions.

    ! Arguments: Input:
    !  grid - grid being searched
    !  icen,jcen - tracker parameter center gridpoints
    !  loncen,latcen - tracker parameter centers' geographic locations
    !  calcparm - is each center valid?
    !  iguess, jguess - first guess gridpoint location
    !  longuess,latguess - first guess geographic location

    ! Arguments: Output:
    !  iout,jout - uv guess center location

    implicit none
    type(fv_atmos_type), intent(inout) :: Atm
    integer, intent(in) :: ids,ide,jds,jde,kds,kde
    integer, intent(in) :: ims,ime,jms,jme,kms,kme
    integer, intent(in) :: its,ite,jts,jte,kts,kte

    integer, intent(in) :: icen(maxtp), jcen(maxtp)
    real, intent(in) :: loncen(maxtp), latcen(maxtp)
    logical, intent(in) :: calcparm(maxtp)

    integer, intent(in) :: iguess,jguess
    real, intent(in) :: latguess,longuess

    integer, intent(inout) :: iout,jout
    real :: degrees,dist
    integer :: ip,ict
    integer(kind=8) :: isum,jsum

    ict=2
    isum=2*iguess
    jsum=2*jguess

    ! Get a guess storm center location for searching for the wind centers:
    do ip=1,maxtp
       if ((ip > 2 .and. ip < 7) .or. ip == 10) then
          cycle   ! because 3-6 are for 850 & 700 u & v and 10 is
          ! for surface wind magnitude.
       elseif(calcparm(ip)) then
          call calcdist (longuess,latguess,loncen(ip),latcen(ip),dist,degrees)
          if(dist<uverrmax) then
             if(ip==1 .or. ip==2 .or. ip==11) then
                isum=isum+2*icen(ip)
                jsum=jsum+2*jcen(ip)
                ict=ict+2
             else
                isum=isum+icen(ip)
                jsum=jsum+jcen(ip)
                ict=ict+1
             endif
          endif
       endif
    enddo

    iout=nint(real(isum)/real(ict))
    jout=nint(real(jsum)/real(ict))
  end subroutine get_uv_guess

  subroutine get_uv_center(Atm,orig, &
       iout,jout,rout,calcparm,lonout,latout, &
       dxdymean,cparm, &
       ids,ide,jds,jde,kds,kde, &
       ims,ime,jms,jme,kms,kme, &
       ips,ipe,jps,jpe,kps,kpe, &
       iuvguess,juvguess)

    implicit none

    integer, intent(in) :: iuvguess,juvguess
    type(fv_atmos_type), intent(inout) :: Atm
    character*(*), intent(in) :: cparm
    real, intent(in) :: dxdymean
    real, intent(inout) :: rout
    integer, intent(inout) :: iout,jout
    logical, intent(inout) :: calcparm
    real, intent(inout) :: latout,lonout
    real, intent(in) :: orig(ims:ime,jms:jme)

    integer, intent(in) :: IDS,IDE,JDS,JDE,KDS,KDE
    integer, intent(in) :: IMS,IME,JMS,JME,KMS,KME
    integer, intent(in) :: IPS,IPE,JPS,JPE,KPS,KPE

    integer :: icen,jcen, i,j, istart,istop, jstart,jstop, ierr
    real :: rcen, srsq

    ! Restrict the search area.  By default, we search everywhere except the boundary:
    istart=max(ids+1,ips)
    istop=min(ide-2,ipe)
    jstart=max(jds+1,jps)
    jstop=min(jde-2,jpe)

    ! If the guess location is given, then further restrict the search area:
    istart=max(istart,iuvguess-nint(rads_vmag/(2.*dxdymean)))
    istop=min(istop,iuvguess+nint(rads_vmag/(2.*dxdymean)))
    jstart=max(jstart,juvguess-nint(rads_vmag/(2.*dxdymean)))
    jstop=min(jstop,juvguess+nint(rads_vmag/(2.*dxdymean)))

    srsq=rads_vmag*rads_vmag*1e6

    icen=-99
    jcen=-99
    rcen=9e9
    do j=jstart,jstop
       do i=istart,istop
          if(orig(i,j)<rcen .and. Atm%distsq(i,j)<srsq) then
             rcen=orig(i,j)
             icen=i
             jcen=j
          endif
       enddo
    enddo

    call mp_reduce_minval(rcen,icen,jcen)

    ! Return result:
    resultif: if(icen==-99 .or. jcen==-99) then
       ! No center found.
       calcparm=.false.
    else
       iout=icen
       jout=jcen
       rout=rcen
       calcparm=.true.
       call get_lonlat(Atm,iout,jout,lonout,latout,ierr, &
            ids,ide, jds,jde, kds,kde, &
            ims,ime, jms,jme, kms,kme, &
            ips,ipe, jps,jpe, kps,kpe)
       if(ierr/=0) then
          calcparm=.false.
          return
       endif
    endif resultif
  end subroutine get_uv_center

  subroutine find_center(Atm,orig,srsq, &
       iout,jout,rout,calcparm,lonout,latout, &
       dxdymean,cparm, &
       ids,ide,jds,jde,kds,kde, &
       ims,ime,jms,jme,kms,kme, &
       ips,ipe,jps,jpe,kps,kpe, &
       iuvguess,juvguess,north_hemi)
    ! This routine replaces the gettrk_main functions find_maxmin and
    ! get_uv_center.

    ! Note: Currently, the smoothing is not yet implemented.

    ! Finds the minimum or maximum value of the smoothed version
    ! (smooth) of the given field (orig).  If a center cannot be
    ! found, sets calcparm=.false., otherwise places the longitude in
    ! lonout and latitude in latout, gridpoint location in (iout,jout)

    ! Mandatory arguments:

    ! grid - grid to search
    ! orig - field to search
    ! smooth - smoothed version of the field (smoothed via relax4e)
    ! iout,jout - center location
    ! rout - center value (min MSLP, min wind, max or min zeta, etc.)
    ! calcparm - true if a center was found, false otherwise
    ! lonout,latout - geographic location of the center
    ! dxdymean - mean grid spacing of the entire domain
    ! cparm - which type of field: zeta, hgt, wind, slp
    ! srsq - square of the maximum radius from domain center to search
    ! ids, ..., kpe - grid, memory and patch dimensions

    ! Optional arguments:

    ! iuvguess,juvguess - first guess center location to restrict search
    ! to a subset of the grid.
    ! north_hemi - we're in the northern hemisphere: true or false?

    implicit none

    integer, intent(in), optional :: iuvguess,juvguess
    type(fv_atmos_type), intent(inout) :: Atm
    character*(*), intent(in) :: cparm
    real, intent(in) :: dxdymean, srsq
    real, intent(inout) :: rout
    integer, intent(inout) :: iout,jout
    logical, intent(inout) :: calcparm
    real, intent(inout) :: latout,lonout
    real, intent(in) :: orig(ims:ime,jms:jme)
    real :: smooth(ims:ime,jms:jme)
    character*255 :: message
    logical, optional :: north_hemi

    integer, intent(in) :: IDS,IDE,JDS,JDE,KDS,KDE
    integer, intent(in) :: IMS,IME,JMS,JME,KMS,KME
    integer, intent(in) :: IPS,IPE,JPS,JPE,KPS,KPE

    integer :: icen,jcen,i,j,ismooth,ierr
    real :: rcen, here, sum, mean, cendist, heredist

    integer :: istart,istop,jstart,jstop,itemp

    logical :: findmin

    ! Restrict the search area.  By default, we search everywhere except the boundary:
    istart=max(ids+1,ips)
    istop=min(ide-2,ipe)
    jstart=max(jds+1,jps)
    jstop=min(jde-2,jpe)

    ! If the guess location is given, then further restrict the search area:
    if(present(iuvguess)) then
       istart=max(istart,iuvguess-nint(rads_vmag/(2.*dxdymean)))
       istop=min(istop,iuvguess+nint(rads_vmag/(2.*dxdymean)))
    endif
    if(present(juvguess)) then
       jstart=max(jstart,juvguess-nint(rads_vmag/(2.*dxdymean)))
       jstop=min(jstop,juvguess+nint(rads_vmag/(2.*dxdymean)))
    endif

    ! Currently smooth has not implemented yet
    do j=jps,min(jde-1,jpe)
       do i=ips,min(ide-1,ipe)
          smooth(i,j)=orig(i,j)
       enddo
    enddo

    ! Figure out whether we're finding a min or max:
    ifmin: if(trim(cparm)=='zeta') then
       if(.not.present(north_hemi)) then
          call mpp_error(FATAL, 'When calling module_tracker find_center for zeta, you must specify the hemisphere parameter.')
       endif
       findmin=.not.north_hemi
    elseif(trim(cparm)=='hgt') then
       findmin=.true.
    elseif(trim(cparm)=='slp') then
       findmin=.true.
    elseif(trim(cparm)=='wind') then
       findmin=.true.
    else
100    format('Invalid parameter cparm="',A,'" in module_tracker find_center')
       write(message,100) trim(cparm)
       call mpp_error(FATAL, message)
    endif ifmin

    ! Find the extremum:
    icen=-99
    jcen=-99
    findminmax: if(findmin) then ! Find a minimum
       rcen=9e9
       do j=jstart,jstop
          do i=istart,istop
             if(smooth(i,j)<rcen .and. Atm%distsq(i,j)<srsq) then
                rcen=smooth(i,j)
                icen=i
                jcen=j
             endif
          enddo
       enddo
       call mp_reduce_minval(rcen,icen,jcen)
    else ! Find a maximum
3014   format(A,' maxval i=',I0,' j=',I0,' r=',F0.3)
       rcen=-9e9
       do j=jstart,jstop
          do i=istart,istop
             if(smooth(i,j)>rcen .and. Atm%distsq(i,j)<srsq) then
                rcen=smooth(i,j)
                icen=i
                jcen=j
             endif
          enddo
       enddo
       call mp_reduce_maxval(rcen,icen,jcen)
    endif findminmax

    ! Return result:
    resultif: if(icen==-99 .or. jcen==-99) then
       ! No center found.
       calcparm=.false.
    else
       iout=icen
       jout=jcen
       rout=rcen
       calcparm=.true.
       call get_lonlat(Atm,iout,jout,lonout,latout,ierr, &
            ids,ide, jds,jde, kds,kde, &
            ims,ime, jms,jme, kms,kme, &
            ips,ipe, jps,jpe, kps,kpe)
       if(ierr/=0) then
          calcparm=.false.
          return
       endif
    endif resultif
  end subroutine find_center

  subroutine get_distsq(Atm, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         its,ite,jts,jte,kts,kte)
    ! This computes approximate distances in km from the domain
    ! center of the various points in the domain. It uses the same
    ! computation as used for distsq: the calculation is done in
    ! gridpoint space, approximating the domain as flat.
    ! Point-to-point distances come from Atm%gridstruct%dxa and Atm%gridstruct%dya.
    ! This routine also determines the distance from the tracker
    ! center location to the nearest point in the domain edge.
    implicit none
    type(fv_atmos_type), intent(inout) :: Atm
    character*255 message
    integer, intent(in) :: ids,ide,jds,jde,kds,kde
    integer, intent(in) :: ims,ime,jms,jme,kms,kme
    integer, intent(in) :: its,ite,jts,jte,kts,kte
    integer i,j,cx,cy,ierr
    integer wilbur,harvey ! filler variables for a function call
    real xfar,yfar,far,xshift,max_edge_distsq,clatr,clonr
    real ylat1,ylat2,xlon1,xlon2,mindistsq

    cx=(ide-ids)/2+ids
    cy=(jde-jds)/2+jds

    call get_lonlat(Atm,cx,cy,clonr,clatr,ierr, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         its,ite,jts,jte,kts,kte)
    if(ierr/=0) then
       call mpp_error(FATAL, 'Domain center location is not inside domain.')
    end if

    do j=jts,min(jte,jde-1)
       do i=its,min(ite,ide-1)
          xfar=(i-cx)*Atm%gridstruct%dxa(i,j)
          yfar=(j-cy)*Atm%gridstruct%dya(i,j)
          far = xfar*xfar + yfar*yfar
          Atm%distsq(i,j)=far
       enddo
    enddo

  end subroutine get_distsq

  subroutine get_tracker_distsq(Atm, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         its,ite,jts,jte,kts,kte)
    ! This computes approximate distances in km from the tracker
    ! center of the various points in the domain.  It uses the same
    ! computation as used for distsq: the calculation is done in
    ! gridpoint space, approximating the domain as flat.
    ! Point-to-point distances come from Atm%gridstruct%dxa and Atm%gridstruct%dya.
    ! This routine also determines the distance from the tracker
    ! center location to the nearest point in the domain edge.
    implicit none
    type(fv_atmos_type), intent(inout) :: Atm
    character*255 message
    integer, intent(in) :: ids,ide,jds,jde,kds,kde
    integer, intent(in) :: ims,ime,jms,jme,kms,kme
    integer, intent(in) :: its,ite,jts,jte,kts,kte
    integer i,j,cx,cy,ierr
    integer wilbur,harvey ! filler variables for a function call
    real xfar,yfar,far,xshift,max_edge_distsq,clatr,clonr
    real ylat1,ylat2,xlon1,xlon2,mindistsq

    cx=Atm%tracker_ifix
    cy=Atm%tracker_jfix

    call get_lonlat(Atm,cx,cy,clonr,clatr,ierr, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         its,ite,jts,jte,kts,kte)
    if(ierr/=0) then
       call mpp_error(FATAL, 'Tracker fix location is not inside domain.')
    end if

    do j=jts,min(jte,jde-1)
       do i=its,min(ite,ide-1)
          xfar=(i-cx)*Atm%gridstruct%dxa(i,j)
          yfar=(j-cy)*Atm%gridstruct%dya(i,j)
          far = xfar*xfar + yfar*yfar
          Atm%tracker_distsq(i,j)=far
       enddo
    enddo

    ! Determine angle.  Note that this is mathematical angle, not
    ! compass angle, and is in geographic lat/lon, not rotated
    ! lat/lon.  (Geographic East=0, geographic North=pi/2.)
    xlon1=clonr ; ylat1=clatr
    call clean_lon_lat(xlon1,ylat1)
    xlon1=xlon1*deg_to_rad
    ylat1=ylat1*deg_to_rad
    do j=jts,min(jte,jde-1)
       do i=its,min(ite,ide-1)
          xlon2=Atm%gridstruct%agrid(i,j,1)*rad_to_deg
          ylat2=Atm%gridstruct%agrid(i,j,2)*rad_to_deg
          call clean_lon_lat(xlon2,ylat2)
          xlon2=xlon2*deg_to_rad
          ylat2=ylat2*deg_to_rad
          Atm%tracker_angle(i,j)=atan2(xlon2-xlon1,ylat2-ylat1)
       enddo
    enddo

    ! Determine the distance between the center location and the
    ! domain edge.
    mindistsq=9e19
    if(jts==jds) then
       mindistsq=min(mindistsq,minval(Atm%tracker_distsq(its:min(ite,ide-1),jds)))
    endif
    if(jte==jde) then
       mindistsq=min(mindistsq,minval(Atm%tracker_distsq(its:min(ite,ide-1),jde-1)))
    endif
    if(its==ids) then
       mindistsq=min(mindistsq,minval(Atm%tracker_distsq(ids,jts:min(jte,jde-1))))
    endif
    if(ite==ide) then
       mindistsq=min(mindistsq,minval(Atm%tracker_distsq(ide-1,jts:min(jte,jde-1))))
    endif
    wilbur=1
    harvey=2
    call mp_reduce_minval(mindistsq,wilbur,harvey)

    Atm%tracker_edge_dist=sqrt(mindistsq)

17  format('Min distance from lon=',F9.3,', lat=',F9.3,' to center is ',F19.3)
    !print *,'Min distance from edge to center is ',Atm%tracker_edge_dist
    write(message,17) clonr, clatr, Atm%tracker_edge_dist
    call mpp_error(NOTE, message)
  end subroutine get_tracker_distsq

  subroutine calcdist(rlonb,rlatb,rlonc,rlatc,xdist,degrees)
    ! Copied from gettrk_main.f
    !
    !     ABSTRACT: This subroutine computes the distance between two
    !               lat/lon points by using spherical coordinates to
    !               calculate the great circle distance between the points.
    !                       Figure out the angle (a) between pt.B and pt.C,
    !             N. Pole   then figure out how much of a % of a great
    !               x       circle distance that angle represents.
    !              / \
    !            b/   \     cos(a) = (cos b)(cos c) + (sin b)(sin c)(cos A)
    !            /     \                                             .
    !        pt./<--A-->\c     NOTE: The latitude arguments passed to the
    !        B /         \           subr are the actual lat vals, but in
    !                     \          the calculation we use 90-lat.
    !               a      \                                      .
    !                       \pt.  NOTE: You may get strange results if you:
    !                         C    (1) use positive values for SH lats AND
    !                              you try computing distances across the
    !                              equator, or (2) use lon values of 0 to
    !                              -180 for WH lons AND you try computing
    !                              distances across the 180E meridian.
    !
    !     NOTE: In the diagram above, (a) is the angle between pt. B and
    !     pt. C (with pt. x as the vertex), and (A) is the difference in
    !     longitude (in degrees, absolute value) between pt. B and pt. C.
    !
    !     !!! NOTE !!! -- THE PARAMETER ecircum IS DEFINED (AS OF THE
    !     ORIGINAL WRITING OF THIS SYSTEM) IN KM, NOT M, SO BE AWARE THAT
    !     THE DISTANCE RETURNED FROM THIS SUBROUTINE IS ALSO IN KM.
    !
    implicit none

    real, intent(inout) :: degrees
    real, intent(out) :: xdist
    real, intent(in) :: rlonb,rlatb,rlonc,rlatc
    real, parameter :: dtr = 0.0174532925199433
    real :: distlatb,distlatc,pole,difflon,cosanga,circ_fract
    !
    if (rlatb < 0.0 .or. rlatc < 0.0) then
       pole = -90.
    else
       pole = 90.
    endif
    !
    distlatb = (pole - rlatb) * dtr
    distlatc = (pole - rlatc) * dtr
    difflon  = abs( (rlonb - rlonc)*dtr )
    !
    cosanga = ( cos(distlatb) * cos(distlatc) + &
         sin(distlatb) * sin(distlatc) * cos(difflon))

    !     This next check of cosanga is needed since I have had ACOS crash
    !     when calculating the distance between 2 identical points (should
    !     = 0), but the input for ACOS was just slightly over 1
    !     (e.g., 1.00000000007), due to (I'm guessing) rounding errors.

    if (cosanga > 1.0) then
       cosanga = 1.0
    endif

    degrees    = acos(cosanga) / dtr
    circ_fract = degrees / 360.
    xdist      = circ_fract * ecircum
    !
    !     NOTE: whether this subroutine returns the value of the distance
    !           in km or m depends on the scale of the parameter ecircum.
    !           At the original writing of this subroutine (7/97), ecircum
    !           was given in km.
    !
    return
  end subroutine calcdist

  subroutine get_lonlat(Atm,iguess,jguess,longuess,latguess,ierr, &
       ids,ide, jds,jde, kds,kde, &
       ims,ime, jms,jme, kms,kme, &
       ips,ipe, jps,jpe, kps,kpe)
    ! Returns the latitude (latguess) and longitude (longuess) of the
    ! specified location (iguess,jguess) in the specified grid.
    implicit none
    integer, intent(in) :: &
         ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         ips,ipe, jps,jpe, kps,kpe
    integer, intent(out) :: ierr
    type(fv_atmos_type), intent(inout) :: Atm
    integer, intent(in) :: iguess,jguess
    real, intent(inout) :: longuess,latguess
    real :: weight,zjunk
    integer :: itemp,jtemp

    ierr=0
    zjunk=1
    if(iguess>=ips .and. iguess<=ipe .and. jguess>=jps .and. jguess<=jpe) then
       weight=1
       longuess=Atm%gridstruct%agrid(iguess,jguess,1)*rad_to_deg
       latguess=Atm%gridstruct%agrid(iguess,jguess,2)*rad_to_deg
       itemp=iguess
       jtemp=jguess
    else
       weight=0
       longuess=-999.9
       latguess=-999.9
       itemp=-99
       jtemp=-99
    endif

    call mp_reduce_maxloc(weight,latguess,longuess,zjunk,itemp,jtemp)

    if(itemp==-99 .and. jtemp==-99) then
       ierr=95
    endif
  end subroutine get_lonlat

  subroutine clean_lon_lat(xlon1,ylat1)
    real, intent(inout) :: xlon1,ylat1
    ! This modifies a (lat,lon) pair so that the longitude fits
    ! between [-180,180] and the latitude between [-90,90], taking
    ! into account spherical geometry.
    ! NOTE: inputs and outputs are in degrees
    xlon1=(mod(xlon1+3600.+180.,360.)-180.)
    ylat1=(mod(ylat1+3600.+180.,360.)-180.)
    if(ylat1>90.) then
       ylat1=180.-ylat1
       xlon1=mod(xlon1+360.,360.)-180.
    elseif(ylat1<-90.) then
       ylat1=-180. - ylat1
       xlon1=mod(xlon1+360.,360.)-180.
    endif
  end subroutine clean_lon_lat

  !----------------------------------------------------------------------------------
  ! These two simple routines return an N, S, E or W for the
  ! hemisphere of a latitude or longitude.
  character(1) function get_lat_ns(lat)
    ! This could be written simply as merge('N','S',lat>=0) if F95 allowed
    implicit none
    real :: lat
    if(lat>=0) then
       get_lat_ns='N'
    else
       get_lat_ns='S'
    endif
  end function get_lat_ns
  character(1) function get_lon_ew(lon)
    ! This could be written simply as merge('E','W',lon>=0) if F95 allowed
    implicit none
    real :: lon
    if(lon>=0) then
       get_lon_ew='E'
    else
       get_lon_ew='W'
    endif
  end function get_lon_ew

  subroutine fv_tracker_post_move(Atm)
    ! This updates the tracker i/j fix location and square of the
    ! distance to the tracker center after a nest move.
    type(fv_atmos_type), intent(inout) :: Atm
    integer :: ierr, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe

    ! Get the grid bounds:
    CALL get_ijk_from_domain(Atm,         &
         ids, ide, jds, jde, kds, kde,    &
         ims, ime, jms, jme, kms, kme,    &
         ips, ipe, jps, jpe, kps, kpe    )

    ! Get the i/j center location from the fix location:
    ierr=0
    call get_nearest_lonlat(Atm,Atm%tracker_ifix,Atm%tracker_jfix, &
         ierr,Atm%tracker_fixlon,Atm%tracker_fixlat, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe)

    ! Get the square of the approximate distance to the tracker center
    ! at all points:
    if(ierr==0) &
         call get_tracker_distsq(Atm, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         ips,ipe,jps,jpe,kps,kpe)
  end subroutine fv_tracker_post_move

  subroutine fv_tracker_move(Atm)
    ! This is for
    implicit none
    type(fv_atmos_type), intent(inout) :: Atm
    integer :: cx,cy
    real xdiff,ydiff
    character*255 message

    Atm%mvnest=.false.
    Atm%move_cd_x=0
    Atm%move_cd_y=0

    if(Atm%grid_number==1) then ! Do nothing for MOAD
       return
    endif

    if(Atm%tracker_gave_up) then
       call mpp_error(NOTE,'Not moving: tracker decided the storm dissapated')
       return
    endif

    if(.not.Atm%tracker_havefix) then
       call mpp_error(NOTE,'Not moving: tracker did not find a storm')
       return
    endif

    ! Calcuate domain center indexes
    cx=(Atm%npx-1)/2+1
    cy=(Atm%npy-1)/2+1

    ! Calculate distance in parent grid index space between storm center and
    ! domain center
    xdiff=(Atm%tracker_ifix-real(cx))/Atm%neststruct%refinement
    ydiff=(Atm%tracker_jfix-real(cy))/Atm%neststruct%refinement

    if(xdiff .ge. 1.0) then
       Atm%move_cd_x=1
    else if(xdiff .le. -1.0) then
       Atm%move_cd_x=-1
    else
       Atm%move_cd_x=0
    endif

    if(ydiff .ge. 1.0) then
       Atm%move_cd_y=1
    else if(ydiff .le. -1.0) then
       Atm%move_cd_y=-1
    else
       Atm%move_cd_y=0
    endif

    if(abs(Atm%move_cd_x)>0 .or. abs(Atm%move_cd_y)>0) then
       Atm%mvnest=.true.
       call mpp_error(NOTE,'Moving: tracker center shifted from nest center')
    else
       Atm%mvnest=.false.
       call mpp_error(NOTE,'Not moving: tracker center is near nest center')
    endif
  end subroutine fv_tracker_move

#endif

end module fv_tracker_mod
