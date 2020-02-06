
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

!>@brief The module 'fv_nggps_diags' computes output diagnostics entirely
!! on 3D pressure levels
!>@details The module is designed for applications that process the full
!!3D fields through the NCEP post-processor.

module fv_nggps_diags_mod

! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>kappa, grav, rdgas</td>
!   </tr>
!   <tr>
!     <td>diag_manager_mod</td>
!     <td>register_diag_field, send_data</td>
!   </tr>
!   <tr>
!     <td>field_manager_mod</td>
!     <td>MODEL_ATMOS</td>
!   </tr>
!   <tr>
!     <td>fms_io_mod</td>
!     <td>set_domain, nullify_domain</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_atmos_type</td>
!   </tr>
!   <tr>
!     <td>fv_diagnostics_mod</td>
!     <td>range_check</td>
!   </tr>
!   <tr>
!     <td>mpp_mod</td>
!     <td>mpp_pe, mpp_root_pe</td>
!   </tr>
!   <tr>
!     <td>tracer_manager_mod</td>
!     <td>get_tracer_names, get_number_tracers, get_tracer_index</td>
!   </tr>
! </table>

 use mpp_mod,            only: mpp_pe, mpp_root_pe,FATAL,mpp_error
 use constants_mod,      only: grav, rdgas
 use fms_io_mod,         only: set_domain, nullify_domain
 use time_manager_mod,   only: time_type, get_time
 use diag_manager_mod,   only: register_diag_field, send_data
 use diag_axis_mod,      only: get_axis_global_length, get_diag_axis, get_diag_axis_name
 use diag_data_mod,      only: output_fields, max_output_fields
 use diag_util_mod,      only: find_input_field
 use tracer_manager_mod, only: get_tracer_names, get_number_tracers, get_tracer_index
 use field_manager_mod,  only: MODEL_ATMOS
 use fv_diagnostics_mod, only: range_check, dbzcalc,max_vv,get_vorticity, &
                               max_uh,max_vorticity,bunkers_vector,       &
                               helicity_relative_CAPS,max_vorticity_hy1 
 use fv_arrays_mod,      only: fv_atmos_type
 use mpp_domains_mod,    only: domain1d, domainUG
#ifdef MULTI_GASES
 use multi_gases_mod,  only:  virq
#endif

 implicit none
 private

 real, parameter:: missing_value = -1.e10
 real, parameter:: stndrd_atmos_ps = 101325.
 real, parameter:: stndrd_atmos_lapse = 0.0065

 logical master
 integer :: id_ua, id_va, id_pt, id_delp, id_pfhy, id_pfnh
 integer :: id_w, id_delz, id_diss, id_ps, id_hs, id_dbz, id_omga
 integer :: kstt_ua, kstt_va, kstt_pt, kstt_delp, kstt_pfhy
 integer :: kstt_pfnh, kstt_w, kstt_delz, kstt_diss, kstt_ps,kstt_hs
 integer :: kend_ua, kend_va, kend_pt, kend_delp, kend_pfhy
 integer :: kend_pfnh, kend_w, kend_delz, kend_diss, kend_ps,kend_hs
 integer :: kstt_dbz, kend_dbz, kstt_omga, kend_omga
 integer :: kstt_windvect, kend_windvect
 integer :: id_wmaxup,id_wmaxdn,kstt_wup, kend_wup,kstt_wdn,kend_wdn
 integer :: id_uhmax03,id_uhmin03,id_uhmax25,id_uhmin25,id_maxvort01
 integer :: id_maxvorthy1,kstt_maxvorthy1,kstt_maxvort01,id_ustm
 integer :: kend_maxvorthy1,kend_maxvort01,id_vstm,id_srh01,id_srh03
 integer :: kstt_uhmax03,kstt_uhmin03,kend_uhmax03,kend_uhmin03
 integer :: kstt_uhmax25,kstt_uhmin25,kend_uhmax25,kend_uhmin25
 integer :: kstt_ustm,kstt_vstm,kend_ustm,kend_vstm,kstt_srh01
 integer :: kstt_srh03,kend_srh01,kend_srh03
 integer :: id_maxvort02,kstt_maxvort02,kend_maxvort02
 integer :: isco, ieco, jsco, jeco, npzo, ncnsto
 integer :: isdo, iedo, jsdo, jedo
 integer :: nlevs
 logical :: hydrostatico
 integer, allocatable :: id_tracer(:), all_axes(:)
 integer, allocatable :: kstt_tracer(:), kend_tracer(:)
 real,    allocatable :: ak(:), bk(:)
 character(20),allocatable :: axis_name(:),axis_name_vert(:)

 logical :: module_is_initialized=.false.
 logical :: use_wrtgridcomp_output=.false.
 integer :: sphum, liq_wat, ice_wat       !< GFDL physics
 integer :: rainwat, snowwat, graupel
 real :: vrange(2) = (/ -330.,  330. /)  !< winds
 real :: wrange(2) = (/ -100.,  100. /)  !< vertical wind
 real :: trange(2) = (/  100.,  350. /)  !< temperature
 real :: skrange(2) = (/ -10000000.0,  10000000.0 /)  !< dissipation estimate for SKEB

! file name
 character(len=64) :: file_name = 'gfs_dyn'

! tracers
 character(len=128)   :: tname
 character(len=256)   :: tlongname, tunits

! wrtout buffer
 real(4), dimension(:,:,:), allocatable, target   :: buffer_dyn
 real(4), dimension(:,:,:,:), allocatable, target :: windvect
 real(4), dimension(:,:), allocatable, target     :: psurf
 real, dimension(:,:), allocatable                :: lon, lat
 real, dimension(:,:),allocatable :: up2,dn2,uhmax03,uhmin03
 real, dimension(:,:),allocatable :: uhmax25,uhmin25,maxvort01
 real, dimension(:,:),allocatable :: maxvorthy1,maxvort02
 public :: fv_nggps_diag_init, fv_nggps_diag, fv_nggps_tavg
#ifdef use_WRTCOMP
 public :: fv_dyn_bundle_setup
#endif

contains

 subroutine fv_nggps_diag_init(Atm, axes, Time)
    type(fv_atmos_type), intent(inout), target :: Atm(:)
    integer, intent(in)         :: axes(4)
    type(time_type), intent(in) :: Time
    integer :: n, i, j, nz

    n = 1
    ncnsto = Atm(1)%ncnst
    npzo   = Atm(1)%npz
    isco = Atm(n)%bd%isc; ieco = Atm(n)%bd%iec
    jsco = Atm(n)%bd%jsc; jeco = Atm(n)%bd%jec
    isdo = Atm(n)%bd%isd; iedo = Atm(n)%bd%ied
    jsdo = Atm(n)%bd%jsd; jedo = Atm(n)%bd%jed
    hydrostatico = Atm(n)%flagstruct%hydrostatic

    call set_domain(Atm(1)%domain)  ! Set domain so that diag_manager can access tile information

    sphum   = get_tracer_index (MODEL_ATMOS, 'sphum')
    liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
    ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')

    rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
    snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
    graupel = get_tracer_index (MODEL_ATMOS, 'graupel')

!--------------------------------------------------------------
! Register main prognostic fields: ps, (u,v), t, omega (dp/dt)
!--------------------------------------------------------------
    allocate(id_tracer(ncnsto))
    allocate(kstt_tracer(ncnsto), kend_tracer(ncnsto))
    id_tracer(:) = 0
    kstt_tracer(:) = 0
    kend_tracer(:) = 0

    if (Atm(n)%flagstruct%write_3d_diags) then
!-------------------
! A grid winds (lat-lon)
!-------------------
       id_ua = register_diag_field ( trim(file_name), 'ucomp', axes(1:3), Time,        &
            'zonal wind', 'm/sec', missing_value=missing_value, range=vrange )
       if (id_ua>0) then
         kstt_ua = 1; kend_ua = npzo
         nlevs = nlevs + npzo
       endif

       id_va = register_diag_field ( trim(file_name), 'vcomp', axes(1:3), Time,        &
            'meridional wind', 'm/sec', missing_value=missing_value, range=vrange)
       if (id_va>0) then
         kstt_va = nlevs+1; kend_va = nlevs+npzo
         nlevs = nlevs + npzo
       endif

       if(id_ua>0 .and. id_va>0) then
         kstt_windvect = 1;  kend_windvect = npzo
         allocate(windvect(3,isco:ieco,jsco:jeco,npzo))
         windvect = 0.
       endif

       if( Atm(n)%flagstruct%hydrostatic ) then 
          id_pfhy = register_diag_field ( trim(file_name), 'pfhy', axes(1:3), Time,        &
               'hydrostatic pressure', 'pa', missing_value=missing_value )
          if (id_pfhy>0) then
             kstt_pfhy = nlevs+1; kend_pfhy = nlevs+npzo
             nlevs = nlevs + npzo
          endif
       else
          id_pfnh = register_diag_field ( trim(file_name), 'pfnh', axes(1:3), Time,        &
               'non-hydrostatic pressure', 'pa', missing_value=missing_value )
          if (id_pfnh>0) then
             kstt_pfnh = nlevs+1; kend_pfnh = nlevs+npzo
             nlevs = nlevs + npzo
          endif
          id_w = register_diag_field ( trim(file_name), 'w', axes(1:3), Time,        &
               'vertical wind', 'm/sec', missing_value=missing_value, range=wrange )
          if (id_w>0) then
             kstt_w = nlevs+1; kend_w = nlevs+npzo
             nlevs = nlevs + npzo
          endif
          id_delz = register_diag_field ( trim(file_name), 'delz', axes(1:3), Time,        &
               'height thickness', 'm', missing_value=missing_value )
          if (id_delz>0) then
             kstt_delz = nlevs+1; kend_delz = nlevs+npzo
             nlevs = nlevs + npzo
          endif
       endif

       id_omga = register_diag_field ( trim(file_name), 'omga', axes(1:3), Time,        &
            'Vertical pressure velocity', 'pa/sec', missing_value=missing_value )
       if (id_omga>0) then
          kstt_omga = nlevs+1; kend_omga = nlevs+npzo
          nlevs = nlevs + npzo
       endif

       id_pt   = register_diag_field ( trim(file_name), 'temp', axes(1:3), Time,       &
            'temperature', 'K', missing_value=missing_value, range=trange )
       if (id_pt>0) then
          kstt_pt = nlevs+1; kend_pt = nlevs+npzo
          nlevs = nlevs + npzo
       endif

       id_delp = register_diag_field ( trim(file_name), 'delp', axes(1:3), Time,        &
            'pressure thickness', 'pa', missing_value=missing_value )
       if (id_delp>0) then
          kstt_delp = nlevs+1; kend_delp = nlevs+npzo
          nlevs = nlevs + npzo
       endif

       !--- diagnostic output for skeb: dissipation estimate
       id_diss = register_diag_field ( trim(file_name), 'diss_est', axes(1:3), Time,    &
            'dissipation estimate', 'none', missing_value=missing_value, range=skrange )
       if (id_delp>0) then
          kstt_diss = nlevs+1; kend_diss = nlevs+npzo
          nlevs = nlevs + npzo
       endif

!--------------------
! Tracer diagnostics:
!--------------------
       do i=1, ncnsto
         call get_tracer_names ( MODEL_ATMOS, i, tname, tlongname, tunits )
         id_tracer(i) = register_diag_field ( file_name, trim(tname),  &
                                   axes(1:3), Time, trim(tlongname), &
                                   trim(tunits), missing_value=missing_value)
          if (id_tracer(i)>0) then
             kstt_tracer(i) = nlevs+1; kend_tracer(i) = nlevs+npzo
             nlevs = nlevs + npzo
          endif
       enddo
!
       id_ps = register_diag_field ( trim(file_name), 'ps', axes(1:2), Time,    &
           'surface pressure', 'pa',  missing_value=missing_value )
       if( id_ps > 0) then
          kstt_ps = nlevs+1; kend_ps = nlevs+1
          nlevs = nlevs + 1
          allocate(psurf(isco:ieco,jsco:jeco))
       endif
!
       id_hs = register_diag_field ( trim(file_name), 'hs', axes(1:2), Time,    &
           'surface geopotential height', 'gpm',  missing_value=missing_value )
       if( id_hs > 0) then
          kstt_hs = nlevs+1; kend_hs = nlevs+1
          nlevs = nlevs + 1
       endif
!
       id_dbz = register_diag_field ( trim(file_name), 'reflectivity', axes(1:3), Time,    &
           'Stoelinga simulated reflectivity', 'dBz', missing_value=missing_value)
       if( rainwat > 0 .and. id_dbz > 0) then
          kstt_dbz = nlevs+1; kend_dbz = nlevs+npzo
          nlevs = nlevs + npzo
       endif
       id_ustm = register_diag_field ( trim(file_name), 'ustm',axes(1:2), Time,       &
          'u comp of storm motion', 'm/s', missing_value=missing_value )
       if( id_ustm > 0) then
          kstt_ustm = nlevs+1; kend_ustm = nlevs+1
          nlevs = nlevs + 1
       endif
       id_vstm = register_diag_field ( trim(file_name), 'vstm',axes(1:2), Time,       &
          'v comp of storm motion', 'm/s', missing_value=missing_value )
       if( id_vstm > 0) then
          kstt_vstm = nlevs+1; kend_vstm = nlevs+1
          nlevs = nlevs + 1
       endif

       id_srh01 = register_diag_field ( trim(file_name), 'srh01',axes(1:2), Time,       &
          '0-1km storm rel. helicity', 'm/s**2', missing_value=missing_value )
       if( id_srh01 > 0) then
          kstt_srh01 = nlevs+1; kend_srh01 = nlevs+1
          nlevs = nlevs + 1
       endif
       id_srh03 = register_diag_field ( trim(file_name), 'srh03',axes(1:2), Time,       &
          '0-3km storm rel. helicity', 'm/s**2', missing_value=missing_value )
       if( id_srh03 > 0) then
          kstt_srh03 = nlevs+1; kend_srh03 = nlevs+1
          nlevs = nlevs + 1
       endif
       id_maxvort01 = register_diag_field ( trim(file_name), 'maxvort01',axes(1:2), Time,       &
          'Max hourly 0-1km vert vorticity', '1/s', missing_value=missing_value )
       if( id_maxvort01 > 0) then
          allocate ( maxvort01(isco:ieco,jsco:jeco) )
          kstt_maxvort01 = nlevs+1; kend_maxvort01 = nlevs+1
          nlevs = nlevs + 1
       endif
       id_maxvort02 = register_diag_field ( trim(file_name), 'maxvort02',axes(1:2), Time,       &
          'Max hourly 0-2km vert vorticity', '1/s', missing_value=missing_value )
       if( id_maxvort02 > 0) then
          allocate ( maxvort02(isco:ieco,jsco:jeco) )
          kstt_maxvort02 = nlevs+1; kend_maxvort02 = nlevs+1
          nlevs = nlevs + 1
       endif
       id_maxvorthy1 = register_diag_field ( trim(file_name), 'maxvorthy1',axes(1:2), Time,       &
          'Max hourly hybrid lev1 vert. vorticity', '1/s', missing_value=missing_value )
       if( id_maxvorthy1 > 0) then
          allocate ( maxvorthy1(isco:ieco,jsco:jeco) )
          kstt_maxvorthy1 = nlevs+1; kend_maxvorthy1 = nlevs+1
          nlevs = nlevs + 1
       endif
       id_wmaxup = register_diag_field ( trim(file_name), 'wmaxup',axes(1:2), Time,       &
          'Max hourly updraft velocity', 'm/s', missing_value=missing_value )
       if( id_wmaxup > 0) then
          allocate ( up2(isco:ieco,jsco:jeco) )
          kstt_wup = nlevs+1; kend_wup = nlevs+1
          nlevs = nlevs + 1 
       endif
       id_wmaxdn = register_diag_field ( trim(file_name), 'wmaxdn',axes(1:2), Time,      &
          'Max hourly downdraft velocity', 'm/s', missing_value=missing_value )
!       write (0,*)'id_wmaxdn in fv_nggps=',id_wmaxdn
       if( id_wmaxdn > 0) then
          allocate ( dn2(isco:ieco,jsco:jeco) )
          kstt_wdn = nlevs+1; kend_wdn = nlevs+1
          nlevs = nlevs + 1
       endif
       id_uhmax03 = register_diag_field ( trim(file_name), 'uhmax03',axes(1:2), Time,      &
          'Max hourly max 0-3km updraft helicity', 'm/s**2', missing_value=missing_value )
!       write (0,*)'id_uhmax03 in fv_nggps=',id_uhmax03
       if( .not.Atm(n)%flagstruct%hydrostatic .and. id_uhmax03 > 0 ) then
          allocate ( uhmax03(isco:ieco,jsco:jeco) )
          kstt_uhmax03 = nlevs+1; kend_uhmax03 = nlevs+1
          nlevs = nlevs + 1
       endif
!
       id_uhmin03 = register_diag_field ( trim(file_name), 'uhmin03',axes(1:2), Time,      &
          'Max hourly min 0-3km updraft helicity', 'm/s**2', missing_value=missing_value )
       if( .not.Atm(n)%flagstruct%hydrostatic .and. id_uhmin03 > 0 ) then
          allocate ( uhmin03(isco:ieco,jsco:jeco) )
          kstt_uhmin03 = nlevs+1; kend_uhmin03 = nlevs+1
          nlevs = nlevs + 1
       endif
!
       id_uhmax25 = register_diag_field ( trim(file_name), 'uhmax25',axes(1:2), Time,      &
          'Max hourly max 2-5km updraft helicity', 'm/s**2', missing_value=missing_value )
       if( .not.Atm(n)%flagstruct%hydrostatic .and. id_uhmax25 > 0 ) then
          allocate ( uhmax25(isco:ieco,jsco:jeco) )
          kstt_uhmax25 = nlevs+1; kend_uhmax25 = nlevs+1
          nlevs = nlevs + 1
       endif
!
       id_uhmin25 = register_diag_field ( trim(file_name), 'uhmin25',axes(1:2), Time,      &
          'Max hourly min 2-5km updraft helicity', 'm/s**2', missing_value=missing_value )
       if( .not.Atm(n)%flagstruct%hydrostatic .and. id_uhmin25 > 0 ) then
          allocate ( uhmin25(isco:ieco,jsco:jeco) )
          kstt_uhmin25 = nlevs+1; kend_uhmin25 = nlevs+1
          nlevs = nlevs + 1
       endif
!
       nz = size(atm(1)%ak)
       allocate(ak(nz))
       allocate(bk(nz))
       do i=1,nz
         ak(i) = atm(1)%ak(i)
         bk(i) = atm(1)%bk(i)
       enddo
!      print *,'in ngpps diag init, ak=',ak(1:5),' bk=',bk(1:5)

! get lon,lat information
       if(.not.allocated(lon)) then
         allocate(lon(isco:ieco,jsco:jeco))
         do j=jsco,jeco
           do i=isco,ieco
             lon(i,j) = Atm(n)%gridstruct%agrid(i,j,1)
           enddo
         enddo
!        if(mpp_pe()==mpp_root_pe())print *,'in fv_dyn bundle,lon=',lon(isco,jsco),lon(ieco-2:ieco,jeco-2:jeco)*180./3.14157
       endif
       if(.not.allocated(lat)) then
         allocate(lat(isco:ieco,jsco:jeco))
         do j=jsco,jeco
           do i=isco,ieco
             lat(i,j) = Atm(n)%gridstruct%agrid(i,j,2)
           enddo
         enddo
!        if(mpp_pe()==mpp_root_pe())print *,'in fv_dyn bundle,lat=',lat(isco,jsco),lat(ieco-2:ieco,jeco-2:jeco)*180./3.14157
       endif
    endif
!
!------------------------------------
! use wrte grid component for output
!------------------------------------
      use_wrtgridcomp_output = .false.

 end subroutine fv_nggps_diag_init


 subroutine fv_nggps_diag(Atm, zvir, Time)

    type(fv_atmos_type), intent(inout) :: Atm(:)
    real,                intent(in):: zvir
    type(time_type),     intent(in) :: Time

    integer :: i, j, k, n, ngc, nq, itrac
    logical :: bad_range
    real    :: ptop, allmax
    real, allocatable :: wk(:,:,:), wk2(:,:,:)
    real, dimension(:,:),allocatable :: ustm,vstm,srh01,srh03

    n = 1
    ngc = Atm(n)%ng
    ptop = Atm(n)%ak(1)
    allmax = -20.
    nq = size (Atm(n)%q,4)
    allocate ( wk(isco:ieco,jsco:jeco,npzo) )
    allocate ( wk2(isco:ieco,jsco:jeco,npzo) )
    allocate ( ustm(isco:ieco,jsco:jeco) )
    allocate ( vstm(isco:ieco,jsco:jeco) )
    allocate ( srh01(isco:ieco,jsco:jeco) )
    allocate ( srh03(isco:ieco,jsco:jeco) )
    if ( Atm(n)%flagstruct%range_warn ) then
         call range_check('DELP', Atm(n)%delp, isco, ieco, jsco, jeco, ngc, npzo, Atm(n)%gridstruct%agrid,    &
                           0.01*ptop, 200.E2, bad_range)
         call range_check('UA', Atm(n)%ua, isco, ieco, jsco, jeco, ngc, npzo, Atm(n)%gridstruct%agrid,   &
                           -250., 250., bad_range)
         call range_check('VA', Atm(n)%va, isco, ieco, jsco, jeco, ngc, npzo, Atm(n)%gridstruct%agrid,   &
                           -250., 250., bad_range)
         call range_check('TA', Atm(n)%pt, isco, ieco, jsco, jeco, ngc, npzo, Atm(n)%gridstruct%agrid,   &
                           150., 350., bad_range) !DCMIP ICs have very low temperatures
    endif
    !--- A-GRID WINDS
    if ( .not. allocated(buffer_dyn)) allocate(buffer_dyn(isco:ieco,jsco:jeco,nlevs))
    if(id_ua > 0) call store_data(id_ua, Atm(n)%ua(isco:ieco,jsco:jeco,:), Time, kstt_ua, kend_ua)
    
    if(id_va > 0) call store_data(id_va, Atm(n)%va(isco:ieco,jsco:jeco,:), Time, kstt_va, kend_va)

    !--- set up 3D wind vector
    if(id_ua>0 .and. id_va>0) then
      do k=1,npzo
        do j=jsco,jeco
          do i=isco,ieco
            windvect(1,i,j,k) = Atm(n)%ua(i,j,k)*cos(lon(i,j)) - Atm(n)%va(i,j,k)*sin(lat(i,j))*sin(lon(i,j))
            windvect(2,i,j,k) = Atm(n)%ua(i,j,k)*sin(lon(i,j)) + Atm(n)%va(i,j,k)*sin(lat(i,j))*cos(lon(i,j))
            windvect(3,i,j,k) =                                  Atm(n)%va(i,j,k)*cos(lat(i,j))
          enddo
        enddo
      enddo
    endif

    !--- W (non-hydrostatic)
    if ( .not.Atm(n)%flagstruct%hydrostatic .and. id_w>0  ) then
       call store_data(id_w, Atm(n)%w(isco:ieco,jsco:jeco,:), Time, kstt_w, kend_w)
    endif

    !--- OMGA (non-hydrostatic)
    if ( id_omga>0  ) then
       call store_data(id_omga, Atm(n)%omga(isco:ieco,jsco:jeco,:), Time, kstt_omga, kend_omga)
    endif

    !--- TEMPERATURE
    if(id_pt   > 0) call store_data(id_pt, Atm(n)%pt(isco:ieco,jsco:jeco,:), Time, kstt_pt, kend_pt)

    !--- TRACERS
    do itrac=1, ncnsto
      call get_tracer_names (MODEL_ATMOS, itrac, tname)
      if (id_tracer(itrac) > 0 .and. itrac.gt.nq) then
        call store_data (id_tracer(itrac), Atm(n)%qdiag(isco:ieco,jsco:jeco,:,itrac), Time,  &
                         kstt_tracer(itrac),kend_tracer(itrac) )
      else
        call store_data (id_tracer(itrac), Atm(n)%q(isco:ieco,jsco:jeco,:,itrac), Time,      &
                         kstt_tracer(itrac),kend_tracer(itrac) )
      endif
    enddo

    !--- DELZ (non-hydrostatic)
    if((.not. Atm(n)%flagstruct%hydrostatic) .and. id_delz > 0) then
       do k=1,npzo
         do j=jsco,jeco
           do i=isco,ieco
             wk(i,j,k) = Atm(n)%delz(i,j,k)
           enddo
         enddo
       enddo
       call store_data(id_delz, wk, Time, kstt_delz, kend_delz)
    endif

    !--- PRESSURE (hydrostatic)
    if( Atm(n)%flagstruct%hydrostatic .and. id_pfhy > 0 ) then
       do k=1,npzo
         do j=jsco,jeco
           do i=isco,ieco         
             wk(i,j,k) = 0.5 *(Atm(n)%pe(i,k,j)+Atm(n)%pe(i,k+1,j))
           enddo
         enddo
       enddo
       call store_data(id_pfhy, wk, Time, kstt_pfhy, kend_pfhy)
    endif

#ifdef GFS_PHYS
    !--- DELP
    if(id_delp > 0 .or. ((.not. Atm(n)%flagstruct%hydrostatic) .and. id_pfnh > 0)) then
       do k=1,npzo
         do j=jsco,jeco
           do i=isco,ieco         
             wk(i,j,k) = Atm(n)%delp(i,j,k)*(1.-sum(Atm(n)%q(i,j,k,2:Atm(n)%flagstruct%nwat)))
           enddo
         enddo
       enddo
       call store_data(id_delp, wk, Time, kstt_delp, kend_delp)
    endif

    !--- Surface Pressure (PS)
    ! Re-compute pressure (dry_mass + water_vapor) surface pressure
    if(id_ps > 0) then
      do k=1,npzo
        do j=jsco,jeco
          do i=isco,ieco
            wk(i,j,k) = Atm(n)%delp(i,j,k)*(1.-sum(Atm(n)%q(i,j,k,2:Atm(n)%flagstruct%nwat)))
          enddo
        enddo
      enddo
      do j=jsco,jeco
        do i=isco,ieco
           psurf(i,j) = ptop
           do k=npzo,1,-1
             psurf(i,j)  = psurf(i,j) + wk(i,j,k)
           enddo
        enddo
      enddo
    endif

    !--- PRESSURE (non-hydrostatic)
    if( (.not. Atm(n)%flagstruct%hydrostatic) .and. id_pfnh > 0) then
       do k=1,npzo
         do j=jsco,jeco
           do i=isco,ieco
             wk(i,j,k) = -wk(i,j,k)/(Atm(n)%delz(i,j,k)*grav)*rdgas*Atm(n)%pt(i,j,k)
#ifdef MULTI_GASES
             wk(i,j,k) = wk(i,j,k) * virq(Atm(n)%q(i,j,k,:))
#else
             wk(i,j,k) = wk(i,j,k) * (1.+zvir*Atm(n)%q(i,j,k,sphum))
#endif
           enddo
         enddo
       enddo
       call store_data(id_pfnh, wk, Time, kstt_pfnh, kend_pfnh)
    endif
#else
    !--- DELP
    if(id_delp > 0) call store_data(id_delp, Atm(n)%delp(isco:ieco,jsco:jeco,:), Time, kstt_delp)

    !--- Surface Pressure (PS)
    if( id_ps > 0) then
      do j=jsco,jeco
        do i=isco,ieco
          psurf(i,j) = Atm(n)%ps(i,j)
        enddo
      enddo
    endif

    !--- PRESSURE (non-hydrostatic)
    if( (.not. Atm(n)%flagstruct%hydrostatic) .and. id_pfnh > 0) then
       do k=1,npzo
         do j=jsco,jeco
           do i=isco,ieco
             wk(i,j,k) = -Atm(n)%delp(i,j,k)/(Atm(n)%delz(i,j,k)*grav)*rdgas*Atm(n)%pt(i,j,k)
#ifdef MULTI_GASES
             wk(i,j,k) = wk(i,j,k)*virq(Atm(n)%q(i,j,k,:)
#else
             wk(i,j,k) = wk(i,j,k)*(1.+zvir*Atm(n)%q(i,j,k,sphum))
#endif
           enddo
         enddo
       enddo
       call store_data(id_pfnh, wk, Time, kstt_pfnh, kend_pfnh)
    endif
#endif

    !--- DISS_EST (skeb: dissipation estimate)
    if(id_diss > 0) call store_data(id_diss, Atm(n)%diss_est(isco:ieco,jsco:jeco,:), Time, kstt_diss, kend_diss)
!
    if(id_ps > 0) then
      if(  use_wrtgridcomp_output ) then
        do j=jsco,jeco
          do i=isco,ieco
            wk(i,j,1) = (psurf(i,j)/stndrd_atmos_ps)**(rdgas/grav*stndrd_atmos_lapse)
          enddo
        enddo
      else
        do j=jsco,jeco
          do i=isco,ieco
            wk(i,j,1) = psurf(i,j)
          enddo
        enddo
      endif
!      print *,'in comput ps, i=',isco,'j=',jsco,'psurf=',psurf(isco,jsco),'stndrd_atmos_ps=',stndrd_atmos_ps, &
!       'rdgas=',rdgas,'grav=',grav,'stndrd_atmos_lapse=',stndrd_atmos_lapse,rdgas/grav*stndrd_atmos_lapse
      call store_data(id_ps, wk, Time, kstt_ps, kend_ps)
    endif
    
    if( id_hs > 0) then
      do j=jsco,jeco
        do i=isco,ieco
          wk(i,j,1) = Atm(n)%phis(i,j)/grav
        enddo
      enddo
      call store_data(id_hs, wk, Time, kstt_hs, kend_hs)
    endif

    !--- 3-D Reflectivity field
    if ( rainwat > 0 .and. id_dbz>0) then
      call dbzcalc(Atm(n)%q, Atm(n)%pt, Atm(n)%delp, Atm(n)%peln, Atm(n)%delz, &
                   wk, wk2, allmax, Atm(n)%bd, npzo, Atm(n)%ncnst, Atm(n)%flagstruct%hydrostatic, &
                   zvir, .false., .false., .false., .true. ) ! GFDL MP has constant N_0 intercept
      call store_data(id_dbz, wk, Time, kstt_dbz, kend_dbz)
    endif

    deallocate ( wk )
    !---u and v comp of storm motion, 0-1, 0-3km SRH 
    if ( id_ustm > 0 .or. id_vstm > 0 .or. id_srh01 > 0 .or. id_srh03 > 0) then
      if ( id_ustm > 0 .and. id_vstm > 0 .and. id_srh01 > 0 .and. id_srh03 > 0) then
        call bunkers_vector(isco,ieco,jsco,jeco,ngc,npzo,zvir,sphum,ustm,vstm,    &
                          Atm(n)%ua,Atm(n)%va, Atm(n)%delz, Atm(n)%q,           &
                          Atm(n)%flagstruct%hydrostatic, Atm(n)%pt, Atm(n)%peln,&
                          Atm(n)%phis, grav)

        call helicity_relative_CAPS(isco,ieco,jsco,jeco,ngc,npzo,zvir,sphum,srh01, &
                                 ustm, vstm,Atm(n)%ua, Atm(n)%va, Atm(n)%delz,  &
                                 Atm(n)%q,Atm(n)%flagstruct%hydrostatic,        &
                                 Atm(n)%pt, Atm(n)%peln, Atm(n)%phis, grav, 0., 1.e3)

        call helicity_relative_CAPS(isco,ieco,jsco,jeco,ngc,npzo,zvir,sphum,srh03, &
                                 ustm, vstm,Atm(n)%ua, Atm(n)%va, Atm(n)%delz,  &
                                 Atm(n)%q,Atm(n)%flagstruct%hydrostatic,        &
                                 Atm(n)%pt, Atm(n)%peln, Atm(n)%phis, grav, 0., 3.e3)

        call store_data(id_ustm, ustm, Time, kstt_ustm, kend_ustm)
        call store_data(id_vstm, vstm, Time, kstt_vstm, kend_vstm)
        call store_data(id_srh01, srh01, Time, kstt_srh01, kend_srh01)
        call store_data(id_srh03, srh03, Time, kstt_srh03, kend_srh03)
       else
         print *,'Missing fields in diag_table'
         print *,'Make sure the following are listed in the diag_table under gfs_dyn:'
         print *,'ustm,vstm,srh01,shr03'
         call mpp_error(FATAL, 'Missing fields in diag_table')
         stop
       endif
      endif
        deallocate ( ustm )
        deallocate ( vstm )
        deallocate ( srh01 )
        deallocate ( srh03 )

    !--- max hourly 0-1km vert. vorticity
    if ( id_maxvort01 > 0) then
      call store_data(id_maxvort01, maxvort01, Time, kstt_maxvort01, kend_maxvort01)
    endif
    !--- max hourly 0-2km vert. vorticity
    if ( id_maxvort02 > 0) then
      call store_data(id_maxvort02, maxvort02, Time, kstt_maxvort02, kend_maxvort02)
    endif
    !--- max hourly hybrid lev 1 vert. vorticity 
    if ( id_maxvorthy1 > 0) then
      call store_data(id_maxvorthy1, maxvorthy1, Time, kstt_maxvorthy1, kend_maxvorthy1)
    endif
!   
    !--- max hourly updraft velocity 
    if ( id_wmaxup > 0) then
      call store_data(id_wmaxup, up2, Time, kstt_wup, kend_wup)
    endif
    !--- max hourly downdraft velocity
    if ( id_wmaxdn > 0) then
      call store_data(id_wmaxdn, dn2, Time, kstt_wdn, kend_wdn)
    endif
    !--- max hourly max 0-3km updraft helicity 
    if ( .not.Atm(n)%flagstruct%hydrostatic .and. id_uhmax03 > 0) then
      call store_data(id_uhmax03, uhmax03, Time, kstt_uhmax03, kend_uhmax03)
    endif
!
    !--- max hourly min 0-3km updraft helicity
    if ( .not.Atm(n)%flagstruct%hydrostatic .and. id_uhmin03 > 0) then
      call store_data(id_uhmin03, uhmin03, Time, kstt_uhmin03, kend_uhmin03)
    endif
!   
    !--- max hourly max 2-5km updraft helicity
    if ( .not.Atm(n)%flagstruct%hydrostatic .and. id_uhmax25 > 0) then
      call store_data(id_uhmax25, uhmax25, Time, kstt_uhmax25, kend_uhmax25)
    endif
!
    !--- max hourly min 2-5km updraft helicity
    if ( .not.Atm(n)%flagstruct%hydrostatic .and. id_uhmin25 > 0) then
      call store_data(id_uhmin25, uhmin25, Time, kstt_uhmin25, kend_uhmin25)
    endif

    call nullify_domain()

 end subroutine fv_nggps_diag

 subroutine fv_nggps_tavg(Atm, Time_step_atmos,avg_max_length,zvir)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    type(time_type),     intent(in) :: Time_step_atmos
    real,                intent(in):: zvir
    integer :: i, j, k, n, ngc, nq, itrac
    integer seconds, days, nsteps_per_reset
    logical, save :: first_call=.true.
    real, save :: first_time = 0.
    integer, save :: kdtt = 0
    real :: avg_max_length
    real,dimension(:,:,:),allocatable :: vort 
    n = 1
    ngc = Atm(n)%ng
    nq = size (Atm(n)%q,4)
!
!Check if any of the max hourly fields are being requested otherwise skip
!
    if(id_wmaxup > 0 .or. id_wmaxdn > 0 .or. id_uhmax03 > 0 .or. id_uhmin03 > 0 &
          .or. id_uhmax25 > 0 .or. id_uhmin25 > 0 .or. id_maxvort01 > 0 &
          .or. id_maxvorthy1 > 0 .or. id_maxvort02 > 0) then
!Make sure the group of max hrly fields listed below are ALL present otherwise
!abort
!
     if(id_wmaxup > 0 .and. id_wmaxdn > 0 .and. id_uhmax03 > 0 .and. id_uhmin03 > 0 &
          .and. id_uhmax25 > 0 .and. id_uhmin25 > 0 .and. id_maxvort01 > 0  & 
          .and. id_maxvorthy1 > 0 .and. id_maxvort02 > 0) then
       allocate ( vort(isco:ieco,jsco:jeco,npzo) )
       if (first_call) then
         call get_time (Time_step_atmos, seconds,  days)
         first_time=seconds
         first_call=.false.
         kdtt=0
       endif
       nsteps_per_reset = nint(avg_max_length/first_time)
       do j=jsco,jeco
          do i=isco,ieco
             if(mod(kdtt,nsteps_per_reset)==0)then
                up2(i,j) = -999.
                dn2(i,j) = 999.
                maxvorthy1(i,j)= 0.
                maxvort01(i,j)= 0.
                maxvort02(i,j)= 0.
             endif
          enddo
        enddo
         call get_vorticity(isco,ieco,jsco,jeco,isdo,iedo,jsdo,jedo, &
                           npzo,Atm(n)%u,Atm(n)%v,vort,    &
                           Atm(n)%gridstruct%dx,Atm(n)%gridstruct%dy,&
                           Atm(n)%gridstruct%rarea)
         call max_vorticity_hy1(isco,ieco,jsco,jeco,npzo,vort,maxvorthy1)
         call max_vorticity(isco,ieco,jsco,jeco,ngc,npzo,zvir, &
                           sphum,Atm(n)%delz,Atm(n)%q, &
                           Atm(n)%flagstruct%hydrostatic, &
                           Atm(n)%pt,Atm(n)%peln,Atm(n)%phis,grav, &
                           vort,maxvort01,0., 1.e3)
         call max_vorticity(isco,ieco,jsco,jeco,ngc,npzo,zvir, &
                           sphum,Atm(n)%delz,Atm(n)%q, &
                           Atm(n)%flagstruct%hydrostatic, &
                           Atm(n)%pt,Atm(n)%peln,Atm(n)%phis,grav, &
                           vort,maxvort02,0., 2.e3)
         if( .not.Atm(n)%flagstruct%hydrostatic ) then
            call max_vv(isco,ieco,jsco,jeco,npzo,ngc,up2,dn2,Atm(n)%pe,Atm(n)%w)
            do j=jsco,jeco
               do i=isco,ieco
                  if(mod(kdtt,nsteps_per_reset)==0)then
                    uhmax03(i,j)= 0.
                    uhmin03(i,j)= 0.
                    uhmax25(i,j)= 0.
                    uhmin25(i,j)= 0.
                  endif
               enddo
             enddo

             call max_uh(isco,ieco,jsco,jeco,ngc,npzo,zvir, &
                           sphum,uhmax03,uhmin03,Atm(n)%w,vort,Atm(n)%delz, &
                           Atm(n)%q,Atm(n)%flagstruct%hydrostatic, &
                           Atm(n)%pt,Atm(n)%peln,Atm(n)%phis,grav, &
                           0., 3.e3)
             call max_uh(isco,ieco,jsco,jeco,ngc,npzo,zvir, &
                           sphum,uhmax25,uhmin25,Atm(n)%w,vort,Atm(n)%delz, &
                           Atm(n)%q,Atm(n)%flagstruct%hydrostatic, &
                           Atm(n)%pt,Atm(n)%peln,Atm(n)%phis,grav, &
                           2.e3, 5.e3)
         endif
    kdtt=kdtt+1
    deallocate (vort)
    else
       print *,'Missing max/min hourly field in diag_table'
       print *,'Make sure the following are listed in the diag_table under gfs_dyn:'
       print *,'wmaxup,wmaxdn,uhmax03,uhmin03,uhmax25,uhmin25,maxvort01,maxvort02 and maxvorthy1'
       call mpp_error(FATAL, 'Missing max hourly fields in diag_table')
       stop
    endif
   endif
 end subroutine fv_nggps_tavg
!
 subroutine store_data(id, work, Time, nstt, nend)
   integer, intent(in)         :: id
   integer, intent(in)         :: nstt, nend
   real, intent(in)            :: work(isco:ieco,jsco:jeco,nend-nstt+1)
   type(time_type), intent(in) :: Time
!
   integer k,j,i,kb
   logical :: used
!
   if( id > 0 ) then
     if( use_wrtgridcomp_output ) then
       do k=1,nend-nstt+1
         do j= jsco,jeco
           do i=isco,ieco
             kb = k + nstt - 1
             buffer_dyn(i,j,kb) = work(i,j,k)
           enddo
         enddo
       enddo
     else
       used = send_data(id, work, Time)
     endif
   endif

 end subroutine store_data

#ifdef use_WRTCOMP

 subroutine fv_dyn_bundle_setup(axes, dyn_bundle, fcst_grid, quilting, rc)
!
!-------------------------------------------------------------
!*** set esmf bundle for dyn output fields
!------------------------------------------------------------
!
   use esmf
   use diag_data_mod, ONLY:  diag_atttype
!
   integer, intent(in)         :: axes(:)
   type(ESMF_FieldBundle),intent(inout)        :: dyn_bundle
   type(ESMF_Grid),intent(inout)               :: fcst_grid
   logical,intent(in)                          :: quilting
   integer,intent(out)                         :: rc


!*** local variables
   integer i, j, k, n
   integer num_axes, id, axis_length, direction, edges
   integer num_attributes, num_field_dyn, axis_typ
   character(255) :: units, long_name, cart_name,axis_direct,edgesS
   character(128) :: output_name, output_file, output_file1, dynbdl_name, shydrostatic
   integer currdate(6), idx1
   logical l3Dvector
   type(domain1d) :: Domain
   type(domainuG) :: DomainU
   real,dimension(:),allocatable :: axis_data
   type(diag_atttype),dimension(:),allocatable :: attributes
   character(2) axis_id

   type(ESMF_Field)                            :: field
!
!jwtest
!   integer :: fieldcount
!   character(128) :: fld_outfilename
!   character(128),dimension(:),allocatable      :: fieldnamelist
!   type(ESMF_Field),dimension(:),allocatable    :: fieldlist
!
!------------------------------------------------------------

! initialize RC
   rc = ESMF_SUCCESS

!--- use wrte grid component for output
   use_wrtgridcomp_output = quilting

! data type
   if(.not. allocated(buffer_dyn))allocate(buffer_dyn(isco:ieco,jsco:jeco,nlevs))
   buffer_dyn=0.
   num_field_dyn = 0.
!
! set output files
   call ESMF_FieldBundleGet(dyn_bundle, name=dynbdl_name,rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     line=__LINE__, &
     file=__FILE__)) &
     return  ! bail out
   idx1 = index(dynbdl_name,'_bilinear')
   if(idx1 > 0) then
     output_file = dynbdl_name(1:idx1-1)
   else
     output_file = 'dyn'
   endif
!
!------------------------------------------------------------
!*** add attributes to the bundle such as subdomain limtis,
!*** axes, output time, etc
!------------------------------------------------------------
!
!*** add attributes
   num_axes = size(axes)
   allocate(all_axes(num_axes))
   all_axes(1:num_axes) = axes(1:num_axes)
!   if(mpp_pe()==mpp_root_pe())print *,'in fv_dyn bundle,num_axes=',num_axes, 'axes=',axes
!
!*** add global attributes in the field bundle:
   call ESMF_AttributeAdd(dyn_bundle, convention="NetCDF", purpose="FV3", &
     attrList=(/"hydrostatic", &
                "ncnsto     ", &
                "ak         ", &
                "bk         "/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     line=__LINE__, &
     file=__FILE__)) &
     return  ! bail out
   if (hydrostatico ) then
     shydrostatic = 'hydrostatic'
   else
     shydrostatic = 'non-hydrostatic'
   endif
   call ESMF_AttributeSet(dyn_bundle, convention="NetCDF", purpose="FV3", &
     name="hydrostatic", value=trim(shydrostatic), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     line=__LINE__, &
     file=__FILE__)) &
     return  ! bail out
!
   call ESMF_AttributeSet(dyn_bundle, convention="NetCDF", purpose="FV3", &
     name="ncnsto", value=ncnsto, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     line=__LINE__, &
     file=__FILE__)) &
     return  ! bail out
!
   call ESMF_AttributeSet(dyn_bundle, convention="NetCDF", purpose="FV3", &
     name="ak", valueList=ak, rc=rc)
!    if(mpp_pe()==mpp_root_pe())print *,'in fv_dyn bundle,after add ak, rc=',rc
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     line=__LINE__, &
     file=__FILE__)) &
     return  ! bail out
!
   call ESMF_AttributeSet(dyn_bundle, convention="NetCDF", purpose="FV3", &
     name="bk", valueList=bk, rc=rc)
!    if(mpp_pe()==mpp_root_pe())print *,'in fv_dyn bundle,after add bk, rc=',rc
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     line=__LINE__, &
     file=__FILE__)) &
     return  ! bail out
!
!*** get axis names
   allocate(axis_name(num_axes))
   do id = 1,num_axes
     call get_diag_axis_name( axes(id), axis_name(id))
   enddo
   if( num_axes>2 ) then
     allocate(axis_name_vert(num_axes-2))
     do id=3,num_axes
       axis_name_vert(id-2) = axis_name(id)
     enddo
!     if(mpp_pe()==mpp_root_pe())print *,'in fv_dyn bundle,axis_name_vert=',axis_name_vert
     call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
       attrList=(/"vertical_dim_labels"/), rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
     call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
       name="vertical_dim_labels", valueList=axis_name_vert, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
   endif

   do id = 1,num_axes
     axis_length =  get_axis_global_length(axes(id)) 
     allocate(axis_data(axis_length))
     call get_diag_axis( axes(id), axis_name(id), units, long_name, cart_name, &
                         direction, edges, Domain, DomainU, axis_data,         &
                         num_attributes=num_attributes,              &
                         attributes=attributes)
!
     edgesS=''
     do i = 1,num_axes
       if(axes(i) == edges) edgesS=axis_name(i)
     enddo

!     if(mpp_pe()==mpp_root_pe())print *,'in fv_dyn bundle,id=',id,'edges=',edges,rc, &
!       'num_attributes=',num_attributes,'edgesS=',trim(edgesS)
!
! Add vertical dimension Attributes to Grid
     if( id>2 ) then
!      if(mpp_pe()==mpp_root_pe())print *,' in dyn add grid, axis_name=',     &
!         trim(axis_name(id)),'axis_data=',axis_data
!
! Previous definition using variable-length character arrays violates the Fortran standards.
! While this worked with Intel compilers, it caused the model to crash in different places
! with both GNU and PGI. Compilers should throw an error at compile time, but it seems that
! they can't handle the "trim(...) // ..." expressions.
! The Standard (Fortran 2003) way to do this correctly is to tell the array constructor
! how long to make the fixed array of characters:
!
!   call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3", &
!   attrList=(/ character(128) :: trim(axis_name(id)),trim(axis_name(id))//":long_name", &
!             trim(axis_name(id))//":units", trim(axis_name(id))//":cartesian_axis", &
!             trim(axis_name(id))//":positive", trim(axis_name(id))//":edges"/), rc=rc)
!
! However this fails for GNU and PGI, see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=85547
! The easiest and safest way forward is to define the attributes one by one as it is done
! as it is done below in add_field_to_bundle.
!
      ! Add attributes one by one
      call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
        attrList=(/trim(axis_name(id))/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
        attrList=(/trim(axis_name(id))//":long_name"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
        attrList=(/trim(axis_name(id))//":units"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
        attrList=(/trim(axis_name(id))//":cartesian_axis"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
        attrList=(/trim(axis_name(id))//":positive"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if(trim(edgesS)/='') then
        call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
          attrList=(/trim(axis_name(id))//":edges"/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      ! Set attributes
      call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
        name=trim(axis_name(id)), valueList=axis_data, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
        name=trim(axis_name(id))//":long_name", value=trim(long_name), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
        name=trim(axis_name(id))//":units", value=trim(units), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
        name=trim(axis_name(id))//":cartesian_axis", value=trim(cart_name), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if(direction>0) then
          axis_direct="up"
      else
          axis_direct="down"
      endif
      call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
        name=trim(axis_name(id))//":positive", value=trim(axis_direct), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if(trim(edgesS)/='') then
        call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
          name=trim(axis_name(id))//":edges", value=trim(edgesS), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

     endif

     deallocate(axis_data)
   enddo
!
!*** add esmf fields
   if(id_ua > 0) then
     call find_outputname(trim(file_name),'ucomp',output_name)
!     if(mpp_pe()==mpp_root_pe()) print *,'ucomp output name is ',trim(output_name)
     call add_field_to_bundle(trim(output_name),'zonal wind', 'm/sec', "time: point",   &
          axes(1:3), fcst_grid, kstt_ua,kend_ua, dyn_bundle, output_file, &
          range=vrange, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
!
   if(id_va > 0) then
     call find_outputname(trim(file_name),'vcomp',output_name)
     call add_field_to_bundle(trim(output_name),'meridional wind', 'm/sec', "time: point",   &
          axes(1:3), fcst_grid, kstt_va,kend_va, dyn_bundle, output_file, &
          range=vrange,rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
!
!*** create 3D vector from local u/v winds
   if(id_ua > 0 .and. id_va > 0) then
     output_name = "windvector"
     output_file1 = 'none'
     l3Dvector = .true.
     call add_field_to_bundle(trim(output_name),'3D cartisian wind vector', 'm/sec', "time: point",   &
          axes(1:3), fcst_grid, kstt_windvect,kend_windvect, dyn_bundle, output_file1, range=vrange,   &
          l3Dvector=l3Dvector,rcd=rc)
   endif
!
   if ( .not.hydrostatico ) then
     if( id_w>0  ) then
       call find_outputname(trim(file_name),'w',output_name)
       call add_field_to_bundle(trim(output_name),'vertical wind', 'm/sec', "time: point",   &
            axes(1:3), fcst_grid, kstt_w,kend_w, dyn_bundle, output_file, &
            range=wrange, rcd=rc)
       if(rc==0)  num_field_dyn=num_field_dyn+1
     endif
     if( id_pfnh>0  ) then
       call find_outputname(trim(file_name),'pfnh',output_name)
       call add_field_to_bundle(trim(output_name),'non-hydrostatic pressure', 'pa', "time: point",  &
            axes(1:3), fcst_grid, kstt_pfnh,kend_pfnh, dyn_bundle, output_file, rcd=rc)
       if(rc==0)  num_field_dyn=num_field_dyn+1
     endif
     if( id_delz>0  ) then
       call find_outputname(trim(file_name),'delz',output_name)
       call add_field_to_bundle(trim(output_name),'height thickness', 'm', "time: point",   &
            axes(1:3), fcst_grid, kstt_delz,kend_delz, dyn_bundle, output_file, rcd=rc)
       if(rc==0)  num_field_dyn=num_field_dyn+1
     endif
   else
     if( id_pfhy>0  ) then
       call find_outputname(trim(file_name),'pfhy',output_name)
       call add_field_to_bundle(trim(output_name),'hydrostatic pressure', 'pa', "time: point",   &
            axes(1:3), fcst_grid, kstt_pfhy,kend_pfhy, dyn_bundle, output_file, rcd=rc)
       if(rc==0)  num_field_dyn=num_field_dyn+1
     endif
   endif
!
   if( id_omga>0  ) then
     call find_outputname(trim(file_name),'omga',output_name)
     call add_field_to_bundle(trim(output_name),'Vertical pressure velocity', 'pa/sec', "time: point",   &
          axes(1:3), fcst_grid, kstt_omga,kend_omga, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
!
   if(id_pt > 0) then
     call find_outputname(trim(file_name),'temp',output_name)
     call add_field_to_bundle(trim(output_name),'temperature', 'K', "time: point",   &
          axes(1:3), fcst_grid, kstt_pt,kend_pt, dyn_bundle, output_file, &
          range=trange,rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
!
   if( id_delp > 0) then
     call find_outputname(trim(file_name),'delp',output_name)
     call add_field_to_bundle(trim(output_name),'pressure thickness', 'pa', "time: point",   &
          axes(1:3), fcst_grid, kstt_delp,kend_delp, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
!
! tracers
   do i=1, ncnsto
     call get_tracer_names ( MODEL_ATMOS, i, tname, tlongname, tunits )
     if (id_tracer(i)>0) then
       call find_outputname(trim(file_name),trim(tname),output_name)
       call add_field_to_bundle(trim(output_name),trim(tlongname), trim(tunits), "time: point",   &
            axes(1:3), fcst_grid, kstt_tracer(i),kend_tracer(i), dyn_bundle, output_file, rcd=rc)
       if(rc==0)  num_field_dyn=num_field_dyn+1
     endif
!     if(mpp_pe()==mpp_root_pe())print *,'in fv_dyn bundle,add trac,i=',i,'output_name=',trim(output_name),' rc=',rc
   enddo
!
!
   if( id_ps > 0) then
     call find_outputname(trim(file_name),'ps',output_name)
     call add_field_to_bundle(trim(output_name),'surface pressure', 'pa', "time: point",   &
          axes(1:2), fcst_grid, kstt_ps,kend_ps, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
!
   if( id_hs > 0) then
     call find_outputname(trim(file_name),'hs',output_name)
     call add_field_to_bundle(trim(output_name),'surface geopotential height', 'gpm', "time: point",   &
          axes(1:2), fcst_grid, kstt_hs,kend_hs, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
!
   if(id_dbz > 0) then
     call find_outputname(trim(file_name),'reflectivity',output_name)
!     if(mpp_pe()==mpp_root_pe())print *,'reflectivity, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'Stoelinga simulated reflectivity', 'dBz', "time: point",   &
          axes(1:3), fcst_grid, kstt_dbz,kend_dbz, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
   if(id_ustm > 0 .and. id_vstm > 0 .and. id_srh01 > 0 .and. id_srh03 > 0) then
     call find_outputname(trim(file_name),'ustm',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'u comp. of storm motion, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'u comp of storm motion', 'm/s', "time: point",   &
          axes(1:2), fcst_grid, kstt_ustm,kend_ustm, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1

     call find_outputname(trim(file_name),'vstm',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'v comp. of storm motion, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'v comp of storm motion', 'm/s', "time: point",   &
          axes(1:2), fcst_grid, kstt_vstm,kend_vstm, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1

     call find_outputname(trim(file_name),'srh01',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'0-1km srh, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'0-1km srh', 'm/s**2', "time: point",   &
          axes(1:2), fcst_grid, kstt_srh01,kend_srh01, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1

     call find_outputname(trim(file_name),'srh03',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'0-3km srh, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'0-3km srh', 'm/s**2', "time: point",   &
          axes(1:2), fcst_grid, kstt_srh03,kend_srh03, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif


   if(id_maxvort01 > 0) then
     call find_outputname(trim(file_name),'maxvort01',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'max hourly 0-1km vert. vorticity, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'Max hourly 0-1km vert. vorticity', '1/s', "time: point",   &
          axes(1:2), fcst_grid, kstt_maxvort01,kend_maxvort01, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
   if(id_maxvort02 > 0) then
     call find_outputname(trim(file_name),'maxvort02',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'max hourly 0-2km vert. vorticity, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'Max hourly 0-2km vert. vorticity', '1/s', "time: point",   &
          axes(1:2), fcst_grid, kstt_maxvort02,kend_maxvort02, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
   if(id_maxvorthy1 > 0) then
     call find_outputname(trim(file_name),'maxvorthy1',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'max hourly lev 1 vert. vorticity output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'Max hourly lev 1 vert vort.', '1/s', "time: point",   &
          axes(1:2), fcst_grid, kstt_maxvorthy1,kend_maxvorthy1, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
   if(id_wmaxup > 0) then
     call find_outputname(trim(file_name),'wmaxup',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'max hourly updraft vel, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'Max hourly updraft velocity', 'm/s', "time: point",   &
          axes(1:2), fcst_grid, kstt_wup,kend_wup, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
   if(id_wmaxdn > 0) then
     call find_outputname(trim(file_name),'wmaxdn',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'max hourly downdraft vel, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'Max hourly downdraft velocity', 'm/s', "time: point",   &
          axes(1:2), fcst_grid, kstt_wdn,kend_wdn, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
   if( .not.hydrostatico .and. id_uhmax03 > 0 ) then
     call find_outputname(trim(file_name),'uhmax03',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'max hourly 0-3km updraft helicity, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'Max hourly 0-3km updraft helicity', 'm/s**2', "time: point",   &
          axes(1:2), fcst_grid, kstt_uhmax03,kend_uhmax03, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
   if( .not.hydrostatico .and. id_uhmin03 > 0 ) then
     call find_outputname(trim(file_name),'uhmin03',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'max hourly 0-3km updraft helicity, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'Max hourly 0-3km updraft helicity', 'm/s**2', "time: point",   &
          axes(1:2), fcst_grid, kstt_uhmin03,kend_uhmin03, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
   if( .not.hydrostatico .and. id_uhmax25 > 0 ) then
     call find_outputname(trim(file_name),'uhmax25',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'max hourly 2-5km updraft helicity, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'Max hourly 2-5km updraft helicity', 'm/s**2', "time: point",   &
          axes(1:2), fcst_grid, kstt_uhmax25,kend_uhmax25, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif
   if( .not.hydrostatico .and. id_uhmin25 > 0 ) then
     call find_outputname(trim(file_name),'uhmin25',output_name)
     if(mpp_pe()==mpp_root_pe())print *,'max hourly 2-5km updraft helicity, output name=',trim(output_name)
     call add_field_to_bundle(trim(output_name),'Max hourly 2-5km updraft helicity', 'm/s**2', "time: point",   &
          axes(1:2), fcst_grid, kstt_uhmin25,kend_uhmin25, dyn_bundle, output_file, rcd=rc)
     if(rc==0)  num_field_dyn=num_field_dyn+1
   endif

!jwtest:
!   call ESMF_FieldBundleGet(dyn_bundle, fieldCount=fieldCount, rc=rc)
!   print *,'in dyn_bundle_setup, fieldCount=',fieldCount
!   allocate(fieldnamelist(fieldCount),fieldlist(fieldCount))
!   call ESMF_FieldBundleGet(dyn_bundle, fieldlist=fieldlist,fieldnamelist=fieldnamelist, rc=rc)
!   do i=1,fieldCount
!     call ESMF_AttributeGet(fieldlist(i), convention="NetCDF", purpose="FV3", &
!                    name="output_file", value=fld_outfilename, rc=rc)
!     print *,'in dyn bundle setup, i=',i,' fieldname=',trim(fieldnamelist(i)),' out filename=',trim(fld_outfilename)
!   enddo

 end subroutine fv_dyn_bundle_setup

 subroutine add_field_to_bundle(var_name,long_name,units,cell_methods, axes,dyn_grid, &
                                kstt,kend,dyn_bundle,output_file, range, l3Dvector, rcd)
   use esmf
   implicit none

   character(*), intent(in)             :: var_name, long_name, units, cell_methods
   integer, intent(in)                  :: axes(:)
   type(esmf_grid), intent(in)          :: dyn_grid
   integer, intent(in)                  :: kstt,kend
   type(esmf_fieldbundle),intent(inout) :: dyn_bundle
   character(*), intent(in)             :: output_file
   real, intent(in), optional           :: range(2)
   logical, intent(in), optional        :: l3Dvector
   integer, intent(out), optional       :: rcd
!
!*** local variable
   type(ESMF_Field)         :: field
   type(ESMF_DataCopy_Flag) :: copyflag=ESMF_DATACOPY_REFERENCE
   integer rc, i, j, idx
   real(4),dimension(:,:,:,:),pointer :: temp_r4d
   real(4),dimension(:,:,:),  pointer :: temp_r3d
   real(4),dimension(:,:),    pointer :: temp_r2d
   logical, save :: first=.true.
!
!*** create esmf field  
   if( present(l3Dvector) ) then
     temp_r4d => windvect(1:3,isco:ieco,jsco:jeco,kstt:kend)
     call ESMF_LogWrite('create winde vector esmf field', ESMF_LOGMSG_INFO, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!jw      field = ESMF_FieldCreate(dyn_grid, temp_r4d, datacopyflag=ESMF_DATACOPY_VALUE, 
     field = ESMF_FieldCreate(dyn_grid, temp_r4d, datacopyflag=ESMF_DATACOPY_REFERENCE, &
                            gridToFieldMap=(/2,3/), ungriddedLBound=(/1,kstt/), ungriddedUBound=(/3,kend/), &
                            name="windvector", indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
     call ESMF_LogWrite('create winde vector esmf field', ESMF_LOGMSG_INFO, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

     call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
        attrList=(/"output_file"/), rc=rc)
     call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
        name='output_file',value=trim(output_file),rc=rc)

     call ESMF_FieldBundleAdd(dyn_bundle,(/field/), rc=rc)
     if( present(rcd)) rcd=rc
     return
   else if( kend>kstt ) then
     temp_r3d => buffer_dyn(isco:ieco,jsco:jeco,kstt:kend)
     field = ESMF_FieldCreate(dyn_grid, temp_r3d, datacopyflag=copyflag, &
                            name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
   else if(kend==kstt) then
     temp_r2d => buffer_dyn(isco:ieco,jsco:jeco,kstt)
     field = ESMF_FieldCreate(dyn_grid, temp_r2d, datacopyflag=copyflag, &
                            name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
   endif
!
!*** add field attributes
   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
        attrList=(/"long_name"/), rc=rc)
   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
        name='long_name',value=trim(long_name),rc=rc)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
        attrList=(/"units"/), rc=rc)
   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
        name='units',value=trim(units),rc=rc)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
        attrList=(/"missing_value"/), rc=rc)
   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
        name='missing_value',value=missing_value,rc=rc)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
        attrList=(/"_FillValue"/), rc=rc)
   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
        name='_FillValue',value=missing_value,rc=rc)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
        attrList=(/"cell_methods"/), rc=rc)
   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
        name='cell_methods',value=trim(cell_methods),rc=rc)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
        attrList=(/"output_file"/), rc=rc)
   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
        name='output_file',value=trim(output_file),rc=rc)
!
!*** add vertical coord attribute:
   if( size(axes) > 2) then
     do i=3,size(axes)
       idx=0
       do j=1,size(all_axes)
         if (axes(i)==all_axes(j)) then
           idx=j
           exit
         endif
       enddo
       if (idx>0) then
         call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
           attrList=(/"ESMF:ungridded_dim_labels"/), rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, &
           file=__FILE__)) &
           return  ! bail out
         call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
           name="ESMF:ungridded_dim_labels", valueList=(/trim(axis_name(idx))/), rc=rc)
!         if( first ) then
!             print *,'add axis_name to field,',trim(axis_name(idx))
!         endif
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, &
           file=__FILE__)) &
           return  ! bail out
       endif
     enddo
     first=.false.
   endif

!*** add field into bundle
   call ESMF_FieldBundleAdd(dyn_bundle,(/field/), rc=rc)
   if( present(rcd)) rcd=rc
!
 end subroutine add_field_to_bundle
!-------------------------------------------------------------------------------------
 subroutine find_outputname(module_name, field_name, output_name)
   character(*), intent(in)     :: module_name
   character(*), intent(in)     :: field_name
   character(*), intent(out)    :: output_name
!
   integer i,j,in_num, out_num
   integer tile_count
!
   tile_count=1
   in_num = find_input_field(module_name, field_name, tile_count)
!
   output_name=''
   do i=1, max_output_fields
     if(output_fields(i)%input_field == in_num) then
       output_name=output_fields(i)%output_name
       exit
     endif
   enddo
   if(output_name=='') then
     print *,'Error, cant find out put name, field_name=',trim(field_name),'in_num=',in_num
   endif

 end subroutine find_outputname

#endif

end module fv_nggps_diags_mod
