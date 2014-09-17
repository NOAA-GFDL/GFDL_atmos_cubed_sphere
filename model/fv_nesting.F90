module fv_nesting_mod

   use mpp_domains_mod,     only: mpp_update_domains
   use mpp_domains_mod,     only: mpp_global_field
   use field_manager_mod,   only: MODEL_ATMOS
   use tracer_manager_mod,  only: get_tracer_index
   use fv_sg_mod,           only: neg_adj3
   use mpp_domains_mod,     only: mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain
   use mpp_domains_mod,     only: DGRID_NE, mpp_update_domains, domain2D
   use fv_restart_mod,      only: d2a_setup, d2c_setup
   use mpp_mod,             only: mpp_sync_self, mpp_sync, mpp_send, mpp_recv, mpp_error, FATAL
   use boundary_mod,        only: nested_grid_BC_save, nested_grid_BC, update_coarse_grid
   use fv_mp_mod,           only: is, ie, js, je, isd, ied, jsd, jed, isc, iec, jsc, jec
   use fv_arrays_mod,       only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_nest_BC_type_3D
   use fv_arrays_mod,       only: allocate_fv_nest_BC_type, fv_atmos_type, fv_grid_bounds_type
   use fv_grid_utils_mod,   only: ptop_min, g_sum, cubed_to_latlon, f_p
   use init_hydro_mod,      only: p_var
   use constants_mod,       only: grav, pi, radius, hlv, rdgas    ! latent heat of water vapor
   use fv_mapz_mod,         only: compute_total_energy, mappm, E_Flux_nest
   use fv_timing_mod,       only: timing_on, timing_off
!!! DEBUG CODE
   use mpp_mod,             only: mpp_pe
!!! END DEBUG CODE

implicit none
   logical :: RF_initialized = .false.
   logical :: bad_range
   real, allocatable ::  rf(:), rw(:)
   integer :: kmax=1
   !Arrays for global grid total energy, used for grid nesting
   real, allocatable :: te_2d_coarse(:,:)
   real, allocatable :: dp1_coarse(:,:,:)

private
public :: twoway_nesting, setup_nested_grid_BCs

!---- version number -----
   character(len=128) :: version = '$Id: fv_nesting.F90,v 20.0 2013/12/13 23:04:32 fms Exp $'
   character(len=128) :: tagname = '$Name: tikal_201409 $'

contains

 subroutine setup_nested_grid_BCs(npx, npy, npz, cp_air, zvir, ncnst, sphum,     &
                        u, v, w, pt, delp, delz,q,   &
                        uc, vc, pkz, &
                        nested, inline_q, make_nh, ng, &
                        gridstruct, flagstruct, neststruct, &
                        nest_timestep, tracer_nest_timestep, &
                        domain, bd, proc_in)

   
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(IN) :: cp_air
    real, intent(IN) :: zvir

    integer, intent(IN) :: npx, npy, npz
    integer, intent(IN) :: ncnst, sphum, ng
    logical, intent(IN) :: inline_q, make_nh,nested

    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !  W (m/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: delz(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! height thickness (m)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) ! (uc,vc) mostly used as the C grid winds
    real, intent(inout) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)
    real, intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             ! finite-volume mean pk
    integer, intent(INOUT) :: nest_timestep, tracer_nest_timestep

    type(fv_grid_type), intent(INOUT) :: gridstruct
    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_nest_type), intent(INOUT), target :: neststruct
    type(domain2d), intent(INOUT) :: domain

    logical, intent(INOUT), OPTIONAL :: proc_in

    real, allocatable :: g_dat(:,:,:), pt_coarse(:,:,:), pkz_coarse(:,:,:)
!!$    real :: pkz_south(isd:ied,jsd:0,  npz)
!!$    real :: pkz_north(isd:ied,npy:jed,npz)
!!$    real :: pkz_west (isd:0,  jsd:jed,npz)
!!$    real :: pkz_east (npx:ied,jsd:jed,npz)
    integer :: i,j,k,n,p
    !integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer :: istart, iend
    logical :: process

#ifdef NESTNOCOMM
    logical, save:: initBCs = .true.
#endif

   type(fv_nest_BC_type_3d) :: pkz_BC

    !local pointers
    logical, pointer :: child_grids(:)

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

    child_grids => neststruct%child_grids
    
    if (present(proc_in)) then
       process = proc_in
    else
       process = .true.
    endif

    !IF nested, set up nested grid BCs for time-interpolation
    !(actually applying the BCs is done in dyn_core

    nest_timestep = 0
    if (.not. inline_q) tracer_nest_timestep = 0

#ifdef NESTNOCOMM
    if (neststruct%nested .and. initBCs) call set_BCs_t0(ncnst, flagstruct%hydrostatic, neststruct)
    initBCs = .false.
    return
#endif

    if (neststruct%nested .and. (.not. (neststruct%first_step) .or. make_nh) ) &
         call set_BCs_t0(ncnst, flagstruct%hydrostatic, neststruct) 

    !compute uc/vc for nested-grid BCs
    !!! CLEANUP: if we compute uc/vc here we don't need to do on the first call of c_sw, right?
    if (ANY(neststruct%child_grids)) then
       call timing_on('COMM_TOTAL')
       !!! CLEANUP: could we make this a non-blocking operation?
       call mpp_update_domains(u, v, &
            domain, gridtype=DGRID_NE, complete=.true.)
       call timing_off('COMM_TOTAL')
       do k=1,npz
          call d2c_setup(u(isd,jsd,k),  v(isd,jsd,k),   &
               uc(isd,jsd,k), vc(isd,jsd,k), flagstruct%nord>0, &
               isd,ied,jsd,jed, is,ie,js,je, npx,npy, &
               gridstruct%grid_type, gridstruct%nested, &
               gridstruct%se_corner, gridstruct%sw_corner, &
               gridstruct%ne_corner, gridstruct%nw_corner, &
               gridstruct%rsin_u, gridstruct%rsin_v, &
               gridstruct%cosa_s, gridstruct%rsin2 )
       end do       
    endif

!! Nested grid: receive from parent grid
    if (neststruct%nested) then
       call nested_grid_BC_save(neststruct%nest_domain, &
            neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy,  npz, bd, &
            neststruct%delp_BC, proc_in=.true.)
       do n=1,ncnst
          call nested_grid_BC_save(neststruct%nest_domain, &
               neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy,  npz, bd, &
               neststruct%q_BC(n), proc_in=.true.)
       enddo
#ifndef SW_DYNAMICS
       call nested_grid_BC_save(neststruct%nest_domain, &
            neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy,  npz, bd, &
            neststruct%pt_BC, proc_in=.true.)


       allocate(pkz_coarse(isd:ied,jsd:jed,npz)) !Filling with pkz
       pkz_coarse = 0.
       call allocate_fv_nest_BC_type(pkz_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz,ng,0,0,0,.false.)
       call nested_grid_BC_save(neststruct%nest_domain, &
            neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy,  npz, bd, &
            pkz_BC, proc_in=.true.)

       call setup_pt_BC(neststruct%pt_BC, pkz_BC, neststruct%q_BC(sphum), npx, npy, npz, cp_air, zvir, bd)

       deallocate(pkz_coarse)

       if (.not. flagstruct%hydrostatic) then
          call nested_grid_BC_save(neststruct%nest_domain, &
               neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy,  npz, bd, &
               neststruct%w_BC, proc_in=.true.)
          call nested_grid_BC_save(neststruct%nest_domain, &
               neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy,  npz, bd, &
               neststruct%delz_BC, proc_in=.true.)
       endif
#endif
       call nested_grid_BC_save(neststruct%nest_domain, &
            neststruct%ind_u, neststruct%wt_u, 0, 1,  npx,  npy,  npz, bd, &
            neststruct%u_BC, proc_in=.true.)
       call nested_grid_BC_save(neststruct%nest_domain, &
            neststruct%ind_u, neststruct%wt_u, 0, 1,  npx,  npy,  npz, bd, &
            neststruct%vc_BC, proc_in=.true.)
       call nested_grid_BC_save(neststruct%nest_domain, &
            neststruct%ind_v, neststruct%wt_v, 1, 0,  npx,  npy,  npz, bd, &
            neststruct%v_BC, proc_in=.true.)
       call nested_grid_BC_save(neststruct%nest_domain, &
            neststruct%ind_v, neststruct%wt_v, 1, 0,  npx,  npy,  npz, bd, &
            neststruct%uc_BC, proc_in=.true.)
    endif


!! Coarse grid: send to child grids
#ifndef SW_DYNAMICS
    if (ANY(child_grids)) then
       allocate(pkz_coarse(isd:ied,jsd:jed,npz))
       pkz_coarse(isc:iec, jsc:jec, :) = pkz
    endif
#endif

    do p=1,size(child_grids)
       if (child_grids(p)) then
          call nested_grid_BC_save(delp, neststruct%nest_domain_all(p), 0, 0, &
               1, npx-1, 1, npy-1, npz)
          do n=1,ncnst
             call nested_grid_BC_save(q(:,:,:,n), neststruct%nest_domain_all(p), 0, 0, &
                  1, npx-1, 1, npy-1, npz)
          enddo
#ifndef SW_DYNAMICS
          call nested_grid_BC_save(pt, neststruct%nest_domain_all(p), 0, 0, &
               1, npx-1, 1, npy-1, npz)

          !Working with PKZ is more complicated since it is only defined on the interior of the grid.
          call nested_grid_BC(pkz_coarse, neststruct%nest_domain_all(p), 0, 0)

          if (.not. flagstruct%hydrostatic) then
             call nested_grid_BC_save(w, neststruct%nest_domain_all(p), 0, 0, &
                  1, npx-1, 1, npy-1, npz)
             call nested_grid_BC_save(delz, neststruct%nest_domain_all(p), 0, 0, &
                  1, npx-1, 1, npy-1, npz)
          endif
          
#endif
          call nested_grid_BC_save(u, neststruct%nest_domain_all(p), 0, 1, &
               1, npx-1, 1, npy-1, npz)
          call nested_grid_BC_save(vc, neststruct%nest_domain_all(p), 0, 1, &
               1, npx-1, 1, npy-1, npz)
          call nested_grid_BC_save(v, neststruct%nest_domain_all(p), 1, 0, &
               1, npx-1, 1, npy-1, npz)
          call nested_grid_BC_save(uc, neststruct%nest_domain_all(p), 1, 0, &
               1, npx-1, 1, npy-1, npz)
       endif
    enddo
    
#ifndef SW_DYNAMICS
    if (ANY(child_grids)) then
       deallocate(pkz_coarse)
    endif
#endif


    if (neststruct%first_step) then
       if (neststruct%nested) call set_BCs_t0(ncnst, flagstruct%hydrostatic, neststruct)
       neststruct%first_step = .false.
    endif

    call mpp_sync_self

 end subroutine setup_nested_grid_BCs

 subroutine setup_pt_BC(pt_BC, pkz_BC, sphum_BC, npx, npy, npz, cp_air, zvir, bd)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(IN), target    :: pkz_BC, sphum_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pt_BC
   integer, intent(IN) :: npx, npy, npz
   real, intent(IN) :: cp_air, zvir

   real, dimension(:,:,:), pointer :: ptBC, pkzBC, sphumBC

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (is == 1) then
      ptBC    =>    pt_BC%west_t1
      pkzBC   =>   pkz_BC%west_t1
      sphumBC => sphum_BC%west_t1
      do k=1,npz
      do j=jsd,jed
      do i=isd,0
         ptBC(i,j,k) = cp_air*ptBC(i,j,k)/pkzBC(i,j,k)*(1.+zvir*sphumBC(i,j,k))
      end do
      end do
      end do
   end if

   if (js == 1) then
      ptBC    =>    pt_BC%south_t1
      pkzBC   =>   pkz_BC%south_t1
      sphumBC => sphum_BC%south_t1
      if (is == 1) then
         istart = is
      else
         istart = isd
      end if
      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      do k=1,npz
      do j=jsd,0
      do i=istart,iend
         ptBC(i,j,k) = cp_air*ptBC(i,j,k)/pkzBC(i,j,k) * &
              (1.+zvir*sphumBC(i,j,k))
      end do
      end do
      end do
   end if


   if (ie == npx-1) then
      ptBC    =>    pt_BC%east_t1
      pkzBC   =>   pkz_BC%east_t1
      sphumBC => sphum_BC%east_t1
      do k=1,npz
      do j=jsd,jed
      do i=npx,ied
         ptBC(i,j,k) = cp_air*ptBC(i,j,k)/pkzBC(i,j,k) * &
              (1.+zvir*sphumBC(i,j,k))
      end do
      end do
      end do
   end if

   if (je == npy-1) then
      ptBC    =>    pt_BC%north_t1
      pkzBC   =>   pkz_BC%north_t1
      sphumBC => sphum_BC%north_t1
      if (is == 1) then
         istart = is
      else
         istart = isd
      end if
      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      do k=1,npz
      do j=npy,jed
      do i=istart,iend
         ptBC(i,j,k) = cp_air*ptBC(i,j,k)/pkzBC(i,j,k) * &
              (1.+zvir*sphumBC(i,j,k))
      end do
      end do
      end do
   end if
   
 end subroutine setup_pt_BC


!!$ subroutine setup_a_nested_grid_BC(var, parent_var, neststruct%nest_domain, istag, jstag, npx, npy, npz, child_grids, BC)
!!$
!!$   real, dimension(:,:,:) :: var, parent_var
!!$   type(neststruct%nest_domain_type), intent(INOUT) :: neststruct%nest_domain
!!$   integer, intent(IN) :: istag, jstag, npx, npy, npz
!!$   logical, intent(IN) :: child_grids
!!$   type(fv_nest_BC_type_3d) :: BC
!!$
!!$
!!$      !Coarse grid: send to child grids
!!$      do p=1,size(child_grids)
!!$         if (child_grids(p)) then
!!$            call nested_grid_BC_save(var, neststruct%nest_domain_all(p), istag,jstag, &
!!$                 1, npx-1, 1, npy-1, npz)
!!$         endif
!!$      enddo
!!$
!!$      !nested grid: receive from parent grids
!!$      if (nested) then
!!$
!!$         call nested_grid_BC_save(parent_var, neststruct%nest_domain, &
!!$              neststruct%ind_v, neststruct%wt_v, 1, 0,  npx,  npy,  npz, 1, parent_grid%npx-1, 1, parent_grid%npy-1, &
!!$              BC, ns=0, proc_in=.true.)
!!$      endif
!!$
!!$
!!$ end subroutine setup_a_nested_grid_BC
!!$

 subroutine set_BCs_t0(ncnst, hydrostatic, neststruct)

   integer, intent(IN) :: ncnst
   logical, intent(IN) :: hydrostatic
   type(fv_nest_type), intent(INOUT) :: neststruct

   integer :: n

   neststruct%delp_BC%east_t0  = neststruct%delp_BC%east_t1
   neststruct%delp_BC%west_t0  = neststruct%delp_BC%west_t1
   neststruct%delp_BC%north_t0 = neststruct%delp_BC%north_t1
   neststruct%delp_BC%south_t0 = neststruct%delp_BC%south_t1
   do n=1,ncnst
      neststruct%q_BC(n)%east_t0  = neststruct%q_BC(n)%east_t1
      neststruct%q_BC(n)%west_t0  = neststruct%q_BC(n)%west_t1
      neststruct%q_BC(n)%north_t0 = neststruct%q_BC(n)%north_t1
      neststruct%q_BC(n)%south_t0 = neststruct%q_BC(n)%south_t1
   enddo
#ifndef SW_DYNAMICS
   neststruct%pt_BC%east_t0    = neststruct%pt_BC%east_t1
   neststruct%pt_BC%west_t0    = neststruct%pt_BC%west_t1
   neststruct%pt_BC%north_t0   = neststruct%pt_BC%north_t1
   neststruct%pt_BC%south_t0   = neststruct%pt_BC%south_t1
   neststruct%pt_BC%east_t0    = neststruct%pt_BC%east_t1
   neststruct%pt_BC%west_t0    = neststruct%pt_BC%west_t1
   neststruct%pt_BC%north_t0   = neststruct%pt_BC%north_t1
   neststruct%pt_BC%south_t0   = neststruct%pt_BC%south_t1
   if (.not. hydrostatic) then
      neststruct%w_BC%east_t0  = neststruct%w_BC%east_t1
      neststruct%w_BC%west_t0  = neststruct%w_BC%west_t1
      neststruct%w_BC%north_t0 = neststruct%w_BC%north_t1
      neststruct%w_BC%south_t0 = neststruct%w_BC%south_t1

      neststruct%delz_BC%east_t0  = neststruct%delz_BC%east_t1
      neststruct%delz_BC%west_t0  = neststruct%delz_BC%west_t1
      neststruct%delz_BC%north_t0 = neststruct%delz_BC%north_t1
      neststruct%delz_BC%south_t0 = neststruct%delz_BC%south_t1

   endif
#endif
   neststruct%u_BC%east_t0  = neststruct%u_BC%east_t1
   neststruct%u_BC%west_t0  = neststruct%u_BC%west_t1
   neststruct%u_BC%north_t0 = neststruct%u_BC%north_t1
   neststruct%u_BC%south_t0 = neststruct%u_BC%south_t1
   neststruct%v_BC%east_t0  = neststruct%v_BC%east_t1
   neststruct%v_BC%west_t0  = neststruct%v_BC%west_t1
   neststruct%v_BC%north_t0 = neststruct%v_BC%north_t1
   neststruct%v_BC%south_t0 = neststruct%v_BC%south_t1


   neststruct%vc_BC%east_t0  = neststruct%vc_BC%east_t1
   neststruct%vc_BC%west_t0  = neststruct%vc_BC%west_t1
   neststruct%vc_BC%north_t0 = neststruct%vc_BC%north_t1
   neststruct%vc_BC%south_t0 = neststruct%vc_BC%south_t1
   neststruct%uc_BC%east_t0  = neststruct%uc_BC%east_t1
   neststruct%uc_BC%west_t0  = neststruct%uc_BC%west_t1
   neststruct%uc_BC%north_t0 = neststruct%uc_BC%north_t1
   neststruct%uc_BC%south_t0 = neststruct%uc_BC%south_t1


 end subroutine set_BCs_t0


!! nestupdate types
!! 1 - Interpolation update on all variables
!! 2 - Conserving update (over areas on cell-
!!     centered variables, over faces on winds) on all variables
!! 3 - Interpolation update on winds only
!! 4 - Interpolation update on all variables except delp (mass conserving)
!! 5 - Remap interpolating update, delp not updated
!! 6 - Remap conserving update, delp not updated
!! 7 - Remap conserving update, delp and q not updated
!! 8 - Remap conserving update, only winds updated

!! Note that nestupdate > 3 will not update delp.

!! "Remap update" remaps updated variables from the nested grid's
!!  vertical coordinate to that of the coarse grid. When delp is not
!!  updated (nestbctype >= 3) the vertical coordinates differ on
!!  the two grids, because the surface pressure will be different
!!  on the two grids.
!! Note: "conserving updates" do not guarantee global conservation
!!  unless flux nested grid BCs are specified, or if a quantity is
!!  not updated at all. For tracers this requires setting
!!  nestbctype > 1; this is not implemented for air mass (delp)

subroutine twoway_nesting(Atm, ngrids, grids_on_this_pe, kappa, cp_air, zvir, dt_atmos)

   type(fv_atmos_type), intent(INOUT) :: Atm(ngrids)
   integer, intent(IN) :: ngrids
   logical, intent(IN) :: grids_on_this_pe(ngrids)
   real, intent(IN) :: kappa, cp_air, zvir, dt_atmos

   integer :: n, p, sphum

   
   if (ngrids > 1) then
      if (Atm(1)%neststruct%parent_of_twoway .and. grids_on_this_pe(1)) then
         call before_twoway_nest_update(Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, Atm(1)%ng, &
              kappa, cp_air, zvir, Atm(1)%ncnst,   &
              Atm(1)%u, Atm(1)%v, Atm(1)%w, Atm(1)%delz, Atm(1)%pt, Atm(1)%delp, Atm(1)%q,   &
              Atm(1)%ps, Atm(1)%pe, Atm(1)%pk, Atm(1)%peln, Atm(1)%pkz, &
              Atm(1)%phis, Atm(1)%ua, Atm(1)%va, &
              Atm(1)%grid_number, Atm(1)%gridstruct, Atm(1)%flagstruct, Atm(1)%idiag, Atm(1)%bd)
      endif

      do n=ngrids,2,-1 !loop backwards to allow information to propagate from finest to coarsest grids

         !two-way updating    
         if (Atm(n)%neststruct%twowaynest ) then
            if  (grids_on_this_pe(n) .or. grids_on_this_pe(Atm(n)%parent_grid%grid_number)) then
               sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
               call twoway_nest_update(Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, cp_air, zvir, &
                    Atm(n)%ncnst, sphum, Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%omga, &
                    Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%uc, Atm(n)%vc, &
                    kappa, Atm(n)%pkz, Atm(n)%delz, Atm(n)%ps, Atm(n)%ptop, &
                    Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%neststruct, Atm(n)%parent_grid, Atm(N)%bd, .false.)
            endif
         endif

      end do

      !NOTE: these routines need to be used with any grid which has been updated to, not just the coarsest grid.
      do n=1,ngrids
         if (Atm(n)%neststruct%parent_of_twoway .and. grids_on_this_pe(n)) then
            call after_twoway_nest_update( Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%ng, dt_atmos,  &
                 kappa, cp_air, zvir, Atm(n)%ncnst,   &
                 Atm(n)%u,  Atm(n)%v,  Atm(n)%w,  Atm(n)%delz, &
                 Atm(n)%pt,  Atm(n)%delp,  Atm(n)%q,   &
#ifdef PKC
                 Atm(n)%ps,  Atm(n)%pe,  Atm(n)%pk,  Atm(n)%peln,  Atm(n)%pkz, Atm(n)%pkc,  &
#else
                 Atm(n)%ps,  Atm(n)%pe,  Atm(n)%pk,  Atm(n)%peln,  Atm(n)%pkz, &
#endif
                 Atm(n)%phis,  Atm(n)%omga,  Atm(n)%ua,  Atm(n)%va,  Atm(n)%uc,  Atm(n)%vc,          &
                 Atm(n)%ptop, Atm(n)%ak,  Atm(n)%bk, Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%idiag, &
                 Atm(n)%ze0,  Atm(n)%grid_number, Atm(n)%domain, Atm(n)%bd)
         endif
      enddo

   endif ! ngrids > 1




  end subroutine twoway_nesting

!!!CLEANUP: this routine assumes that the PARENT GRID has pt = (regular) temperature,
!!!not potential temperature; which may cause problems when updating if this is not the case.
 subroutine twoway_nest_update(npx, npy, npz, cp_air, zvir, ncnst, sphum,     &
                        u, v, w, omga, pt, delp, q,   &
                        uc, vc, kappa, pkz, delz, ps, ptop, &
                        gridstruct, flagstruct, neststruct, &
                        parent_grid, bd, conv_theta_in)

    real, intent(IN) :: kappa, cp_air
    real, intent(IN) :: zvir, ptop

    integer, intent(IN) :: npx, npy, npz
    integer, intent(IN) :: ncnst, sphum
    logical, intent(IN), OPTIONAL :: conv_theta_in

    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !  W (m/s)
    real, intent(inout) :: omga(bd%isd:bd%ied,bd%jsd:bd%jed,npz)      ! Vertical pressure velocity (pa/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) ! (uc,vc) C grid winds
    real, intent(inout) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)

    real, intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             ! finite-volume mean pk
    real, intent(inout) :: delz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! delta-height (m); non-hydrostatic only
    real, intent(inout) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)           ! Surface pressure (pascal)

    type(fv_grid_type), intent(INOUT) :: gridstruct
    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_nest_type), intent(INOUT) :: neststruct

    type(fv_atmos_type), intent(INOUT) :: parent_grid

    real, allocatable :: g_dat(:,:,:)
    real, allocatable :: pt_coarse(:,:,:), pkz_coarse(:,:,:)
    real, allocatable :: t_nest(:,:,:), ps0(:,:)
    integer :: i,j,k,n, r, s
    integer :: isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p
    integer :: isg, ieg, jsg,jeg, npx_p, npy_p
    integer :: istart, iend
    integer :: upoff = 1 !Want 1 for H-S runs?
    !integer :: upoff = 0 !Want 1 for H-S runs?
    real :: rg
    logical :: used, conv_theta=.true.

    real :: qdp(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real, allocatable :: qdp_coarse(:,:,:)
    real(kind=f_p), allocatable :: q_diff(:,:,:)
    
    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    if (neststruct%nestbctype > 1) upoff = 0

    !We update actual temperature, not theta.
    !If pt is actual temperature, set conv_theta to .false.
    if (present(conv_theta_in)) conv_theta = conv_theta_in

    rg = kappa*cp_air
    
    if ((.not. neststruct%parent_proc) .and. (.not. neststruct%child_proc)) return

    call mpp_get_data_domain( parent_grid%domain, &
         isd_p,  ied_p,  jsd_p,  jed_p  )
    call mpp_get_compute_domain( parent_grid%domain, &
         isc_p,  iec_p,  jsc_p,  jec_p  )

   r = neststruct%refinement
   s = r/2 !rounds down (since r > 0)

   !delp/ps

   if (neststruct%nestupdate < 3) then

         call update_coarse_grid(parent_grid%delp, delp, neststruct%nest_domain,&
              neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
              isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npx, npy, npz, 0, 0, &
              neststruct%refinement, neststruct%nestupdate, upoff, 0, neststruct%parent_proc, neststruct%child_proc, parent_grid)

      call mpp_sync!self

#ifdef SW_DYNAMICS
      if (neststruct%parent_proc) then
         do j=jsd_p,jed_p
            do i=isd_p,ied_p

               parent_grid%ps(i,j) = &
                    parent_grid%delp(i,j,1)/grav 

            end do
         end do
      endif
#endif

   end if

   !if (neststruct%nestupdate /= 3 .and. neststruct%nestbctype /= 3) then
   if (neststruct%nestupdate /= 3 .and. neststruct%nestupdate /= 7 .and. neststruct%nestupdate /= 8) then

      allocate(qdp_coarse(isd_p:ied_p,jsd_p:jed_p,npz))
      if (parent_grid%flagstruct%nwat > 0) then
         allocate(q_diff(isd_p:ied_p,jsd_p:jed_p,npz))
         q_diff = 0.
      endif

      do n=1,ncnst

         if (neststruct%child_proc) then
            do k=1,npz
            do j=jsd,jed
            do i=isd,ied
               qdp(i,j,k) = q(i,j,k,n)*delp(i,j,k)
            enddo
            enddo
            enddo
         else
            qdp = 0.
         endif

         if (neststruct%parent_proc) then
            do k=1,npz
            do j=jsd_p,jed_p
            do i=isd_p,ied_p
               qdp_coarse(i,j,k) = parent_grid%q(i,j,k,n)*parent_grid%delp(i,j,k)
            enddo
            enddo
            enddo
            if (n <= parent_grid%flagstruct%nwat) then
            do k=1,npz
            do j=jsd_p,jed_p
            do i=isd_p,ied_p
               q_diff(i,j,k) = q_diff(i,j,k) - parent_grid%q(i,j,k,n)
               !parent_grid%delp(i,j,k) = parent_grid%delp(i,j,k) - qdp_coarse(i,j,k)
            enddo
            enddo
            enddo
            endif
         else
            qdp_coarse = 0.
         endif

            call update_coarse_grid(qdp_coarse, qdp, neststruct%nest_domain, &
                 neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npx, npy, npz, 0, 0, &
                 neststruct%refinement, neststruct%nestupdate, upoff, 0, neststruct%parent_proc, neststruct%child_proc, parent_grid)

               call mpp_sync!self

         if (neststruct%parent_proc) then
            do k=1,npz
            do j=jsd_p,jed_p
            do i=isd_p,ied_p
               parent_grid%q(i,j,k,n) = qdp_coarse(i,j,k)/parent_grid%delp(i,j,k)
            enddo
            enddo
            enddo
            if (n <= parent_grid%flagstruct%nwat) then
            do k=1,npz
            do j=jsd_p,jed_p
            do i=isd_p,ied_p
               q_diff(i,j,k) = q_diff(i,j,k) + parent_grid%q(i,j,k,n)
               !parent_grid%delp(i,j,k) = parent_grid%delp(i,j,k) + qdp_coarse(i,j,k)
            enddo
            enddo
            enddo
            endif
         endif

      end do

         if (neststruct%parent_proc) then
            if (parent_grid%flagstruct%nwat > 0) then
               do k=1,npz
!!$               do j=jsc_p+1,jec_p-1
!!$               do i=isc_p+1,iec_p-1
               do j=jsd_p,jed_p
               do i=isd_p,ied_p
                  parent_grid%delp(i,j,k) = parent_grid%delp(i,j,k)*(1. + q_diff(i,j,k))
               enddo
               enddo
               enddo
               do n=1,parent_grid%flagstruct%nwat
               do k=1,npz
               do j=jsd_p,jed_p
               do i=isd_p,ied_p
!!$               do j=jsc_p+1,jed_p-1
!!$               do i=isc_p+1,ied_p-1
                  parent_grid%q(i,j,k,n) = parent_grid%q(i,j,k,n)/(1. + q_diff(i,j,k))
               enddo
               enddo
               enddo               
               enddo
            endif
         endif

      deallocate(qdp_coarse)
      if (parent_grid%flagstruct%nwat > 0) deallocate(q_diff)

   endif

#ifndef SW_DYNAMICS
   if (neststruct%nestupdate /= 3 .and. neststruct%nestupdate /= 8) then

      if (conv_theta) then

         if (neststruct%child_proc) then
            !pt is potential temperature on the nested grid, but actual
            !temperature on the coarse grid. Compute actual temperature
            !on the nested grid, then gather.
            allocate(t_nest(isd:ied,jsd:jed,1:npz))
            do k=1,npz
               do j=js,je
                  do i=is,ie
                     t_nest(i,j,k) = pt(i,j,k)*pkz(i,j,k)/(cp_air*(1.+zvir*q(i,j,k,sphum)))
                  enddo
               enddo
            enddo
         endif

            call update_coarse_grid(parent_grid%pt, &
                 t_nest, neststruct%nest_domain, &
                 neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npx, npy, npz, 0, 0, &
                 neststruct%refinement, neststruct%nestupdate, upoff, 0, neststruct%parent_proc, neststruct%child_proc, parent_grid)
      else

            call update_coarse_grid(parent_grid%pt, &
              pt, neststruct%nest_domain, &
              neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
              isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npx, npy, npz, 0, 0, &
              neststruct%refinement, neststruct%nestupdate, upoff, 0, neststruct%parent_proc, neststruct%child_proc, parent_grid)

      endif !conv_theta

      call mpp_sync!self

      if (.not. flagstruct%hydrostatic) then

            call update_coarse_grid(parent_grid%w, w, neststruct%nest_domain, &
                 neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npx, npy, npz, 0, 0, &
                 neststruct%refinement, neststruct%nestupdate, upoff, 0, neststruct%parent_proc, neststruct%child_proc, parent_grid)
            !Updating for delz not yet implemented; may be problematic
!!$            call update_coarse_grid(parent_grid%delz, delz, neststruct%nest_domain, &
!!$                 neststruct%ind_update_h, &
!!$                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npz, 0, 0, &
!!$                 neststruct%refinement, neststruct%nestupdate, upoff, 0, neststruct%parent_proc, neststruct%child_proc)

         call mpp_sync!self

      end if
      
   end if !Neststruct%nestupdate /= 3

#endif

      call update_coarse_grid(parent_grid%u, u, neststruct%nest_domain, &
           neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
           isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npx, npy, npz, 0, 1, &
           neststruct%refinement, neststruct%nestupdate, upoff, 0, neststruct%parent_proc, neststruct%child_proc, parent_grid)

      call update_coarse_grid(parent_grid%v, v, neststruct%nest_domain, &
           neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
           isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npx, npy, npz, 1, 0, &
           neststruct%refinement, neststruct%nestupdate, upoff, 0, neststruct%parent_proc, neststruct%child_proc, parent_grid)

   call mpp_sync!self


#ifndef SW_DYNAMICS
   if (neststruct%nestupdate >= 5 .and. npz > 4) then

      !Use PS0 from nested grid, NOT the full delp. Also we assume the same number of levels on both grids.
      !PS0 should be initially set to be ps so that this routine does NOTHING outside of the update region

      !Re-compute nested (AND COARSE) grid ps

      allocate(ps0(isd_p:ied_p,jsd_p:jed_p))
      if (neststruct%parent_proc) then

         parent_grid%ps = parent_grid%ptop
         do k=1,npz
            do j=jsd_p,jed_p
               do i=isd_p,ied_p
                  parent_grid%ps(i,j) = parent_grid%ps(i,j) + &
                       parent_grid%delp(i,j,k)
               end do
            end do
         end do

         ps0 = parent_grid%ps
      endif

      if (neststruct%child_proc) then

         ps = ptop
         do k=1,npz
            do j=jsd,jed
               do i=isd,ied
                  ps(i,j) = ps(i,j) + delp(i,j,k)
               end do
            end do
         end do
      endif

      call update_coarse_grid(ps0, ps, neststruct%nest_domain, &
              neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
              isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npx, npy, 0, 0, &
              neststruct%refinement, neststruct%nestupdate, upoff, 0, neststruct%parent_proc, neststruct%child_proc, parent_grid)

      !!! The mpp version of update_coarse_grid does not return a consistent value of ps
      !!! across PEs, as it does not go into the haloes of a given coarse-grid PE. This
      !!! update_domains call takes care of the problem.

   if (neststruct%parent_proc) then
     call mpp_update_domains(parent_grid%ps, parent_grid%domain, complete=.false.)
     call mpp_update_domains(ps0, parent_grid%domain, complete=.true.)
   endif


      call mpp_sync!self

      if (parent_grid%tile == neststruct%parent_tile) then 

         if (neststruct%parent_proc) then

         !comment out if statement to always remap theta instead of t in the remap-update.
         !(In LtE typically we use remap_t = .true.: remapping t is better (except in
         !idealized simulations with a background uniform theta) since near the top
         !boundary theta is exponential, which is hard to accurately interpolate with a spline
         if (.not. parent_grid%flagstruct%remap_t) then
            do k=1,npz
               do j=jsc_p,jec_p
                  do i=isc_p,iec_p
                     parent_grid%pt(i,j,k) = &
                          cp_air*parent_grid%pt(i,j,k)/parent_grid%pkz(i,j,k)*&
                          (1.+zvir*parent_grid%q(i,j,k,sphum))
                  end do
               end do
            end do
         end if
         !This does not yet do anything with q
         call update_remap_tq(npz, parent_grid%ak, parent_grid%bk, &
              parent_grid%ps, parent_grid%delp, &
              parent_grid%pt, parent_grid%q, npz, ps0, zvir, parent_grid%ptop, ncnst, &
              isc_p, iec_p, jsc_p, jec_p, isd_p, ied_p, jsd_p, jed_p)
         if (.not. parent_grid%flagstruct%remap_t) then
            do k=1,npz
               do j=jsc_p,jec_p
                  do i=isc_p,iec_p
                     parent_grid%pt(i,j,k) = &
                          parent_grid%pt(i,j,k)*parent_grid%pkz(i,j,k) / &
                          (cp_air*(1.+zvir*parent_grid%q(i,j,k,sphum)))
                  end do
               end do
            end do
         end if

         call update_remap_uv(npz, parent_grid%ak, parent_grid%bk, &
              parent_grid%ps, &
              parent_grid%u, &
              parent_grid%v, npz, ps0, &
              isc_p, iec_p, jsc_p, jec_p, isd_p, ied_p, jsd_p, jed_p, parent_grid%ptop)

         endif !neststruct%parent_proc

      end if

      if (allocated(ps0)) deallocate(ps0)

   end if

#endif

 end subroutine twoway_nest_update

!Do ua,va need to be converted back from lat-lon coords?
 subroutine before_twoway_nest_update(npx, npy, npz, ng, &
                        kappa, cp_air, zvir, ncnst,   &
                        u, v, w, delz, pt, delp, q,   &
                        ps, pe, pk, peln, pkz, phis, ua, va, &
                        grid_number, gridstruct, flagstruct, idiag, bd)

    real, intent(IN) :: kappa, cp_air
    real, intent(IN) :: zvir

    integer, intent(IN) :: npx, npy, npz
    integer, intent(IN) :: ng
    integer, intent(IN) :: ncnst

    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !  W (m/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: delz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! delta-height (m); non-hydrostatic only

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)           ! Surface pressure (pascal)
    real, intent(inout) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)  ! edge pressure (pascal)
    real, intent(inout) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)          ! pe**cappa
    real, intent(inout) :: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)           ! ln(pe)
    real, intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             ! finite-volume mean pk
    
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, intent(inout) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)       ! Surface geopotential (g*Z_surf)

    real, intent(inout), dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz):: ua, va

    integer, intent(IN) :: grid_number

    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_grid_type),  intent(INOUT) :: gridstruct
    type(fv_diag_type), intent(IN) :: idiag

    real::   teq(bd%is:bd%ie,bd%js:bd%je)
    integer i, j, k, iq
    real rg

    integer :: sphum, liq_wat, ice_wat
    logical :: used
    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed


    sphum = get_tracer_index(MODEL_ATMOS, 'sphum')
    if (flagstruct%nwat >= 3) then
       liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
       ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
    else
       liq_wat = -1
       ice_wat = -1
    endif

   !Calculate cubed-sphere a-grid winds
   
    !!! CLEANUP: Is this necessary? It is needed for compute_total_energy
    do k=1,npz
       call d2a_setup(u(isd,jsd,k), v(isd,jsd,k), ua(isd,jsd,k), va(isd,jsd,k), .true., &
            isd,ied,jsd,jed, is,ie,js,je, npx,npy, &
            flagstruct%grid_type, gridstruct%nested, gridstruct%cosa_s, gridstruct%rsin2)
    enddo

   if (grid_number /= 1) return   

#ifndef SW_DYNAMICS

   if (.not. allocated(te_2d_coarse)) allocate(te_2d_coarse(isc:iec, jsc:jec))
   if (.not. allocated(dp1_coarse)) allocate(dp1_coarse(is:ie,js:je,npz))
   do k=1,npz
      do j=js,je
         do i=is,ie
            dp1_coarse(i,j,k) = zvir*q(i,j,k,sphum)
         enddo
      enddo
   enddo
   
   rg = kappa*cp_air
   call compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, npz,  &
        u, v, w, delz, pt, delp, q, dp1_coarse, pe, &
        peln, phis, gridstruct%rsin2, gridstruct%cosa_s, zvir, cp_air,  rg, hlv, te_2d_coarse, &
        ua, va, teq, flagstruct%moist_phys, sphum, liq_wat, ice_wat, flagstruct%hydrostatic,idiag%id_te)

#endif

 end subroutine before_twoway_nest_update

 subroutine after_twoway_nest_update(npx, npy, npz, ng, bdt,               &
                        kappa, cp_air, zvir, ncnst,   &
                        u, v, w, delz, pt, delp, q,   &
                        ps, pe, pk, peln, pkz, phis, omga, ua, va, uc, vc,          &
                        ptop, ak, bk, gridstruct, flagstruct, idiag, &
                        ze0, grid_number, domain, bd)

   type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(IN) :: bdt  ! Large time-step
    real, intent(IN) :: kappa, cp_air, ptop
    real, intent(IN) :: zvir

    integer, intent(IN) :: ng, npx, npy, npz
    integer, intent(IN) :: ncnst

    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !  W (m/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: delz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! delta-height (m); non-hydrostatic only
    real, intent(inout) ::  ze0(bd%is:bd%ie,bd%js:bd%je,npz+1) ! height at edges (m); non-hydrostatic

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)           ! Surface pressure (pascal)
    real, intent(inout) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)  ! edge pressure (pascal)
    real, intent(inout) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)          ! pe**cappa
    real, intent(inout) :: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)           ! ln(pe)
    real, intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             ! finite-volume mean pk
    
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, intent(inout) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)       ! Surface geopotential (g*Z_surf)
    real, intent(inout) :: omga(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! Vertical pressure velocity (pa/s)
    real, intent(inout) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) ! (uc,vc) mostly used as the C grid winds
    real, intent(inout) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)

    real, intent(inout), dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz):: ua, va
    real, intent(in),    dimension(npz+1):: ak, bk
    type(fv_grid_type), intent(IN) :: gridstruct
    type(fv_flags_type), intent(IN) :: flagstruct
    type(fv_diag_type), intent(IN) :: idiag
    type(domain2d), intent(INOUT) :: domain

    integer, intent(IN) :: grid_Number

    real :: akap, tpe, rg
    integer:: kord_tracer(ncnst), cld_amt, iq
    real te_2d_coarse_after(bd%is:bd%ie,bd%js:bd%je)
    
    integer :: sphum, liq_wat, ice_wat, i, j, k
    real::   teq(bd%is:bd%ie,bd%js:bd%je)

    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed
    
    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    sphum = get_tracer_index(MODEL_ATMOS, 'sphum')
    if (flagstruct%nwat >= 3) then
       liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
       ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
    else
       liq_wat = -1
       ice_wat = -1
    endif

    cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
    akap  = kappa

    call cubed_to_latlon(u, v, ua, va, &
         gridstruct, npx, npy, npz, &
         1, gridstruct%grid_type, domain, &
         gridstruct%nested, flagstruct%c2l_ord, bd)

#ifndef SW_DYNAMICS

   !To get coarse grid pkz, etc right after a two-way update so
   !that it is consistent across a restart:
   !(should only be called after doing such an update)

    !! CLEANUP: move to twoway_nest_update??
   call p_var(npz, is, ie, js, je, ptop, ptop_min,  &
        delp, delz, &
        pt, ps, &
        pe, peln,   &
        pk,   pkz, kappa, &
        q, ng, flagstruct%ncnst,  gridstruct%area, 0.,  &
        .false.,  .false., & !mountain argument not used
        flagstruct%moist_phys,  flagstruct%hydrostatic, &
        flagstruct%nwat, domain, .false.)

   !DIAGNOSTIC: compute energy after update
   if (grid_number /= 1) return   

   do k=1,npz
      do j=js,je
         do i=is,ie
            dp1_coarse(i,j,k) = zvir*q(i,j,k,sphum)
         enddo
      enddo
   enddo
   
   rg = kappa*cp_air
   call compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, npz,  &
        u, v, w, delz, pt, delp, q, dp1_coarse, pe, &
        peln, phis, gridstruct%rsin2, gridstruct%cosa_s, zvir, cp_air,  rg, hlv, te_2d_coarse_after, &
        ua, va, teq, flagstruct%moist_phys, sphum, liq_wat, ice_wat, flagstruct%hydrostatic,idiag%id_te)

   te_2d_coarse = te_2d_coarse - te_2d_coarse_after
   tpe = g_sum(domain, te_2d_coarse, is, ie, js, je, ng, gridstruct%area, 0)
   E_Flux_nest = tpe / (grav*bdt*4.*pi*radius**2)

#endif

 end subroutine after_twoway_nest_update

 !Routines for remapping (interpolated) nested-grid data to the coarse-grid's vertical coordinate.

 !This does not yet do anything for the tracers
 subroutine update_remap_tq( npz, ak,  bk,  ps, delp,  t,  q,  &
                      kmd, ps0, zvir, ptop, nq, &
                      is, ie, js, je, isd, ied, jsd, jed)
  integer, intent(in):: npz, kmd, nq
  real,    intent(in):: zvir, ptop
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in), dimension(isd:ied,jsd:jed):: ps0
  real,    intent(in), dimension(isd:ied,jsd:jed):: ps
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real,    intent(inout), dimension(isd:ied,jsd:jed,npz):: t
  real,    intent(in), dimension(isd:ied,jsd:jed,npz,nq):: q
  integer,  intent(in) ::  is, ie, js, je, isd, ied, jsd, jed
! local:
  real, dimension(is:ie,kmd):: tp, qp
  real, dimension(is:ie,kmd+1):: pe0, pn0
  real, dimension(is:ie,npz):: qn1
  real, dimension(is:ie,npz+1):: pe1, pn1
  integer i,j,k,n

  do 5000 j=js,je

     do k=1,kmd+1
        do i=is,ie
           pe0(i,k) = ak(k) + bk(k)*ps0(i,j)
           pn0(i,k) = log(pe0(i,k))
       enddo
     enddo 
!------
! Model
!------
     do k=1,kmd+1
        do i=is,ie
           pe1(i,k) = ak(k) + bk(k)*ps(i,j)
           pn1(i,k) = log(pe1(i,k))
       enddo
     enddo 

   do k=1,kmd
      do i=is,ie
         tp(i,k) = t(i,j,k)
      enddo
   enddo
   call mappm(kmd, pn0, tp, npz, pn1, qn1, is,ie, 1, 8, ptop)

        do k=1,npz
           do i=is,ie
              t(i,j,k) = qn1(i,k)
           enddo
        enddo

5000 continue

 end subroutine update_remap_tq

 !remap_uv as-is remaps only a-grid velocities. A new routine has been written to handle staggered grids.
 subroutine update_remap_uv(npz, ak, bk, ps, u, v, kmd, ps0, &
                      is, ie, js, je, isd, ied, jsd, jed, ptop)
  integer, intent(in):: npz
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in):: ps(isd:ied,jsd:jed)
  real,    intent(inout), dimension(isd:ied,jsd:jed+1,npz):: u
  real,    intent(inout), dimension(isd:ied+1,jsd:jed,npz):: v
!
  integer, intent(in):: kmd
  real,    intent(IN) :: ptop
  real,    intent(in):: ps0(isd:ied,jsd:jed)
  integer,  intent(in) ::  is, ie, js, je, isd, ied, jsd, jed
!
! local:
  real, dimension(is:ie+1,kmd+1):: pe0
  real, dimension(is:ie+1,npz+1):: pe1
  real, dimension(is:ie+1,kmd):: qt
  real, dimension(is:ie+1,npz):: qn1
  integer i,j,k

!------
! map u
!------
  do j=js,je+1
!------
! Data
!------
     do k=1,kmd+1
       do i=is,ie
          pe0(i,k) = ak(k) + bk(k)*0.5*(ps0(i,j)+ps0(i,j-1))
       enddo
     enddo
!------
! Model
!------
     do k=1,kmd+1
        do i=is,ie
          pe1(i,k) = ak(k) + bk(k)*0.5*(ps(i,j)+ps(i,j-1))
       enddo
     enddo
!------
!Do map
!------
     qt = 0.
      do k=1,kmd
         do i=is,ie
            qt(i,k) = u(i,j,k)
         enddo
      enddo
      qn1 = 0.
      call mappm(kmd, pe0(is:ie,:), qt(is:ie,:), npz, pe1(is:ie,:), qn1(is:ie,:), is,ie, -1, 8, ptop)
      do k=1,npz
         do i=is,ie
            u(i,j,k) = qn1(i,k)
         enddo
      enddo

   end do

!------
! map v
!------
   do j=js,je

!------
! Data
!------
     do k=1,kmd+1
        do i=is,ie+1
          pe0(i,k) = ak(k) + bk(k)*0.5*(ps0(i,j)+ps0(i-1,j))
       enddo
     enddo
!------
! Model
!------
     do k=1,kmd+1
        do i=is,ie+1
          pe1(i,k) = ak(k) + bk(k)*0.5*(ps(i,j)+ps(i-1,j))
       enddo
     enddo
!------
!Do map
!------
     qt = 0.
      do k=1,kmd
         do i=is,ie+1
            qt(i,k) = v(i,j,k)
         enddo
      enddo
      qn1 = 0.
      call mappm(kmd, pe0(is:ie+1,:), qt(is:ie+1,:), npz, pe1(is:ie+1,:), qn1(is:ie+1,:), is,ie+1, -1, 8, ptop)
      do k=1,npz
         do i=is,ie+1
            v(i,j,k) = qn1(i,k)
         enddo
      enddo
   end do

 end subroutine update_remap_uv


end module fv_nesting_mod
