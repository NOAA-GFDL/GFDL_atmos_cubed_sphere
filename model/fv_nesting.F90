module fv_nesting_mod

   use mpp_domains_mod,     only: mpp_update_domains
   use mpp_domains_mod,     only: mpp_global_field
   use field_manager_mod,   only: MODEL_ATMOS
   use tracer_manager_mod,  only: get_tracer_index
   use fv_sg_mod,           only: neg_adj3
   use mpp_domains_mod,     only: mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain
   use fv_restart_mod,      only: setup_nested_BC_ucvc
   use fv_restart_mod,      only: d2a_setup
   use fv_current_grid_mod, only: grid_type
   use fv_current_grid_mod, only: nest_domain, pelist, child_grids, master, nestbctype
   use fv_current_grid_mod, only: idiag, nested, current_Atm, parent_grid, parent_tile, nest_timestep, tracer_nest_timestep, nsponge, refinement, nestupdate, ng, inline_q
   use mpp_mod,             only: mpp_sync_self, mpp_sync, mpp_send, mpp_recv
   use boundary_mod,        only: nested_grid_BC_save, nested_grid_BC, update_coarse_grid
   use fv_current_grid_mod, only: ind_h, wt_h, ind_u, wt_u, ind_v, wt_v, first_step, ind_update_h
   use fv_current_grid_mod, only: h_east, h_west, h_south, h_north
   use fv_current_grid_mod, only: q_east, q_west, q_south, q_north
#ifndef SW_DYNAMICS
   use fv_current_grid_mod, only: pt_east, pt_west, pt_south, pt_north
   use fv_current_grid_mod, only: w_east, w_west, w_south, w_north
#endif
   use fv_current_grid_mod, only: u_east, u_west, u_south, u_north
   use fv_current_grid_mod, only: v_east, v_west, v_south, v_north
   use fv_current_grid_mod, only: uc_east_t1, uc_west_t1, uc_south_t1, uc_north_t1
   use fv_current_grid_mod, only: vc_east_t1, vc_west_t1, vc_south_t1, vc_north_t1
   use fv_current_grid_mod, only: uc_east_t0, uc_west_t0, uc_south_t0, uc_north_t0
   use fv_current_grid_mod, only: vc_east_t0, vc_west_t0, vc_south_t0, vc_north_t0

   use fv_current_grid_mod, only: h_east_t0, h_west_t0, h_south_t0, h_north_t0
   use fv_current_grid_mod, only: q_east_t0, q_west_t0, q_south_t0, q_north_t0
#ifndef SW_DYNAMICS
   use fv_current_grid_mod, only: pt_east_t0, pt_west_t0, pt_south_t0, pt_north_t0
   use fv_current_grid_mod, only: w_east_t0, w_west_t0, w_south_t0, w_north_t0
#endif
   use fv_current_grid_mod, only: u_east_t0, u_west_t0, u_south_t0, u_north_t0
   use fv_current_grid_mod, only: v_east_t0, v_west_t0, v_south_t0, v_north_t0
   use fv_current_grid_mod, only: uc_east_t0, uc_west_t0, uc_south_t0, uc_north_t0
   use fv_current_grid_mod, only: vc_east_t0, vc_west_t0, vc_south_t0, vc_north_t0
   use fv_current_grid_mod, only: is, ie, js, je, isd, ied, jsd, jed, isc, iec, jsc, jec
   use fv_current_grid_mod, only: hydrostatic, parent_proc, child_proc, rsin_u, rsin_v, area, moist_phys, nwat, cosa_s, rsin2
   use fv_arrays_mod,       only: Atm
   use fv_grid_utils_mod,   only: ptop_min, g_sum, ptop
   use init_hydro_mod,      only: p_var
   use fv_mp_mod,           only: concurrent, gid
   use constants_mod,       only: grav, pi, radius, hlv, rdgas    ! latent heat of water vapor
   use fv_mapz_mod,         only: compute_total_energy, mappm, E_Flux_nest
   use fv_timing_mod,       only: timing_on, timing_off

implicit none
   logical :: RF_initialized = .false.
   logical :: bad_range
   real, allocatable ::  rf(:), rw(:)
   integer :: kmax=1
   !Arrays for global grid total energy, used for grid nesting
   real, allocatable :: te_2d_coarse(:,:)
   real, allocatable :: dp1_coarse(:,:,:)
private
public :: twoway_nest_update,  before_twoway_nest_update, after_twoway_nest_update, setup_nested_grid_BCs

!---- version number -----
   character(len=128) :: version = '$Id: fv_nesting.F90,v 1.1.4.2 2012/05/14 21:23:54 Lucas.Harris Exp $'
   character(len=128) :: tagname = '$Name: siena_201207 $'

contains

 subroutine setup_nested_grid_BCs(npx, npy, npz, cp_air, zvir, ncnst, sphum,     &
                        u, v, w, hydrostatic, pt, delp, q,   &
                        uc, vc, pkz, proc_in)


    real, intent(IN) :: cp_air
    real, intent(IN) :: zvir

    integer, intent(IN) :: npx, npy, npz
    integer, intent(IN) :: ncnst, sphum
    logical, intent(IN) :: hydrostatic

    real, intent(inout), dimension(isd:ied  ,jsd:jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(isd:ied+1,jsd:jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   isd:ied  ,jsd:jed  ,npz)  !  W (m/s)
    real, intent(inout) :: pt(  isd:ied  ,jsd:jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(isd:ied  ,jsd:jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   isd:ied  ,jsd:jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: uc(isd:ied+1,jsd:jed  ,npz) ! (uc,vc) mostly used as the C grid winds
    real, intent(inout) :: vc(isd:ied  ,jsd:jed+1,npz)
    real, intent(inout) :: pkz (is:ie,js:je,npz)             ! finite-volume mean pk

    logical, intent(INOUT), OPTIONAL :: proc_in

    real, allocatable :: g_dat(:,:,:), pt_coarse(:,:,:), pkz_coarse(:,:,:)
    real, allocatable, dimension(:,:,:) :: pkz_east, pkz_west, pkz_south, pkz_north
    integer :: i,j,k,n,p
    !integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer :: isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p
    integer :: isg, ieg, jsg,jeg, npx_p, npy_p
    integer :: istart, iend
    logical :: process

    if (present(proc_in)) then
       process = proc_in
    else
       process = .true.
    endif

      !IF nested, set up nested grid BCs for time-interpolation
      !(actually applying the BCs is done in dyn_core

         nest_timestep = 0
         if (.not. inline_q) tracer_nest_timestep = 0

         if (.not. concurrent) then
            
            call mpp_get_data_domain( parent_grid%domain, &
                 isd_p,  ied_p,  jsd_p,  jed_p  )
            call mpp_get_compute_domain( parent_grid%domain, &
                 isc_p,  iec_p,  jsc_p,  jec_p  )
            call mpp_get_global_domain( parent_grid%domain, &
                 isg, ieg, jsg, jeg, xsize=npx_p, ysize=npy_p)

         endif

         if (nested .and. .not. (first_step .and. concurrent) ) call set_BCs_t0(ncnst) 

         if (concurrent) then

            !Coarse grid: send to child grids
            do p=1,size(child_grids)
               if (child_grids(p)) then
                                           call timing_on('COMM_TOTAL')
                  call nested_grid_BC_save(delp, Atm(p)%nest_domain, 0, 0, &
                       1, npx-1, 1, npy-1, npz)
                                           call timing_off('COMM_TOTAL')
               endif
            enddo

            !nested grid: receive from parent grids
            if (nested) then

                                           call timing_on('COMM_TOTAL')
               call nested_grid_BC_save(parent_grid%delp, nest_domain, &
                    ind_h, wt_h, 0, 0,  npx,  npy,  npz, 1, parent_grid%npx-1, 1, parent_grid%npy-1, &
                    var_east=h_east, &
                    var_west=h_west, &
                    var_north=h_north, &
                    var_south=h_south, ns=0, proc_in=.true.)
                                           call timing_off('COMM_TOTAL')
            endif

         else

                                           call timing_on('COMM_TOTAL')
            call nested_grid_BC_save(parent_grid%delp, nest_domain, &
              ind_h, wt_h, 0, 0, npx, npy, npz, isg, ieg, jsg, jeg, &
              var_east=h_east, &
              var_west=h_west, &
              var_north=h_north, &
              var_south=h_south, ns=0, proc_in=ANY(pelist == gid))
                                           call timing_off('COMM_TOTAL')

         endif

         do n=1,ncnst
            if (concurrent) then

               !Coarse grid: send to child grids
               do p=1,size(child_grids)
                  if (child_grids(p)) then
                                           call timing_on('COMM_TOTAL')
                     call nested_grid_BC_save(q(:,:,:,n), Atm(p)%nest_domain, 0, 0, &
                          1, npx-1, 1, npy-1, npz)
                                           call timing_off('COMM_TOTAL')
                  endif
               enddo

               !nested grid: receive from parent grids
               if (nested) then

                                           call timing_on('COMM_TOTAL')
                  call nested_grid_BC_save(parent_grid%q(:,:,:,n), nest_domain, &
                       ind_h, wt_h, 0, 0,  npx,  npy,  npz, 1, parent_grid%npx-1, 1, parent_grid%npy-1, &
                       var_east=q_east(:,:,:,n), &
                       var_west=q_west(:,:,:,n), &
                       var_north=q_north(:,:,:,n), &
                       var_south=q_south(:,:,:,n), ns=0, proc_in=.true.)
                       !var_south=q_south(:,:,:,n), ns=0, proc_in=.true.)
                                           call timing_off('COMM_TOTAL')

               endif

            else

                                           call timing_on('COMM_TOTAL')
               call nested_grid_BC_save(parent_grid%q(:,:,:,n), nest_domain, &
                 ind_h, wt_h, 0, 0, npx, npy, npz, isg, ieg, jsg, jeg, &
                 var_east=q_east(:,:,:,n), &
                 var_west=q_west(:,:,:,n), &
                 var_north=q_north(:,:,:,n), &
                 var_south=q_south(:,:,:,n), ns=nsponge, proc_in=ANY(pelist == gid)  )
                                           call timing_off('COMM_TOTAL')

            endif

         end do

#ifndef SW_DYNAMICS
         !We will use the virtual potential temperature for the BCs; however, the
         !coarse grid pt is currently the actual temperature. Again, we want to
         !interpolate pkz and compute theta on the nested grid; so we will apply
         !the interpolated pkz to the saved-boundary variables

         if (concurrent) then

            !Coarse grid: send to child grids
            do p=1,size(child_grids)
               if (child_grids(p)) then
                                           call timing_on('COMM_TOTAL')
                  call nested_grid_BC_save(pt, Atm(p)%nest_domain, 0, 0, &
                       1, npx-1, 1, npy-1, npz)
                                           call timing_off('COMM_TOTAL')
               endif
            enddo

            !nested grid: receive from parent grids
            if (nested) then

                                           call timing_on('COMM_TOTAL')
               call nested_grid_BC_save(parent_grid%pt, nest_domain, &
                    ind_h, wt_h, 0, 0,  npx,  npy,  npz, 1, parent_grid%npx-1, 1, parent_grid%npy-1, &
                    var_east=pt_east, &
                    var_west=pt_west, &
                    var_north=pt_north, &
                    var_south=pt_south, ns=0, proc_in=.true.)
                                           call timing_off('COMM_TOTAL')
            endif            

         else

                                           call timing_on('COMM_TOTAL')
            call nested_grid_BC_save(parent_grid%pt, nest_domain, &
!              ind_h, wt_h, 0, 0, npx, npy, npz, isg, ieg, jsg, jeg, &
              ind_h, wt_h, 0, 0, npx, npy, npz, isd_p, ied_p, jsd_p, jed_p, &
              var_east=pt_east, &
              var_west=pt_west, &
              var_north=pt_north, &
              var_south=pt_south, ns=nsponge, proc_in=ANY(pelist == gid)  )
                                           call timing_off('COMM_TOTAL')

         endif



         if (concurrent) then

            if (ANY(child_grids)) then
               allocate(g_dat(isd:ied,jsd:jed,npz)) !Filling with pkz
               g_dat(isc:iec, jsc:jec, :) = pkz
            endif

            !Coarse grid: send to child grids
            do p=1,size(child_grids)
               if (child_grids(p)) then
                                            call timing_on('COMM_TOTAL')
                 call nested_grid_BC(g_dat, Atm(p)%nest_domain, &
                       0, 0,  1, npx-1, 1, npy-1, npz)
                                           call timing_off('COMM_TOTAL')
               endif
            enddo

            if (ANY(child_grids)) then
               deallocate(g_dat)
            endif

            !nested grid: receive from parent grids
            if (nested) then

               allocate(pkz_coarse(isd:ied,jsd:jed,npz)) !Filling with pkz
               pkz_coarse = 0.
               allocate(g_dat(1,1,1))
               g_dat = 0.

!!$               if (is == 1) then
!!$                  allocate(pkz_west(isd:0,jsd:jed,npz))
!!$               else
!!$                  allocate(pkz_west(1,1,1))
!!$               end if
!!$
!!$               if (js == 1) then
!!$                  allocate(pkz_south(max(1,isd):min(npx-1,ied),jsd:0,npz))
!!$               else
!!$                  allocate(pkz_south(1,1,1))
!!$               end if
!!$
!!$               if (ie == npx-1) then
!!$                  allocate(pkz_east(npx:ied,jsd:jed,npz))
!!$               else
!!$                  allocate(pkz_east(1,1,1))
!!$               endif
!!$
!!$               if (je == npy-1) then
!!$                  allocate(pkz_north(max(1,isd):min(npx-1,ied), npy:jed, npz))
!!$               else
!!$                  allocate(pkz_north(1,1,1))
!!$               endif

                                           call timing_on('COMM_TOTAL')
               call nested_grid_BC(pkz_coarse, parent_grid%pkz, nest_domain, &
                    ind_h, wt_h, &
                    0, 0,  npx, npy, npz, isg, ieg, jsg, jeg) 
                                           call timing_off('COMM_TOTAL')
               deallocate(g_dat)

!!$               call nested_grid_BC_save(parent_grid%pkz, nest_domain, &
!!$                    ind_h, wt_h, 0, 0,  npx,  npy,  npz, 1, parent_grid%npx-1, 1, parent_grid%npy-1, &
!!$                    var_east=pkz_east, &
!!$                    var_west=pkz_west, &
!!$                    var_north=pkz_north, &
!!$                    var_south=pkz_south, ns=0, proc_in=.true.)

               call mpp_sync_self

            endif

         else

!! FIXME: nested_grid_BC does not work very well with pkz, since pkz does not have a halo. 
!!  Why doesn't this work? (Does the mpp nesting routine not like variables without haloes?)
!!$         allocate(g_dat(1,1,1))
!!$         allocate(pkz_coarse(isd:ied,jsd:jed,npz)) !Filling with pkz
!!$         call nested_grid_BC(pkz_coarse, parent_grid%pkz, nest_domain, &
!!$                 ind_h, wt_h, &
!!$                 0, 0,  npx, npy, npz, isc_p, iec_p, jsc_p, jec_p)     

            allocate(g_dat(isd_p:ied_p, jsd_p:jed_p,npz))
            g_dat(isc_p:iec_p, jsc_p:jec_p, :) = parent_grid%pkz
            allocate(pkz_coarse(isd:ied,jsd:jed,npz)) !Filling with pkz
            pkz_coarse = 0.
                                           call timing_on('COMM_TOTAL')
            call nested_grid_BC(pkz_coarse, g_dat, nest_domain, &
                 ind_h, wt_h, &
                 0, 0,  npx, npy, npz, isg, ieg, jsg, jeg)     
                                           call timing_off('COMM_TOTAL')

         endif

         if (nested) then


         if (is == 1) then
            do k=1,npz
            do j=jsd,jed
            do i=isd,0
               pt_west(i,j,k) = cp_air*pt_west(i,j,k)/pkz_coarse(i,j,k)*(1.+zvir*q_west(i,j,k,sphum))
            end do
            end do
            end do
         end if

         if (js == 1) then
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
               pt_south(i,j,k) = cp_air*pt_south(i,j,k)/pkz_coarse(i,j,k) * &
                    (1.+zvir*q_south(i,j,k,sphum))
            end do
            end do
            end do
         end if


         if (ie == npx-1) then
            do k=1,npz
            do j=jsd,jed
            do i=npx,ied
               pt_east(i,j,k) = cp_air*pt_east(i,j,k)/pkz_coarse(i,j,k) * &
                    (1.+zvir*q_east(i,j,k,sphum))
            end do
            end do
            end do
         end if

         if (je == npy-1) then
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
               pt_north(i,j,k) = cp_air*pt_north(i,j,k)/pkz_coarse(i,j,k) * &
                    (1.+zvir*q_north(i,j,k,sphum))
            end do
            end do
            end do
         end if

         deallocate(pkz_coarse)

!!$         endif
!!$
         endif

         if (.not. hydrostatic) then

            if (concurrent) then

               !Coarse grid: send to child grids
               do p=1,size(child_grids)
                  if (child_grids(p)) then
                                           call timing_on('COMM_TOTAL')
                     call nested_grid_BC_save(w, Atm(p)%nest_domain, 0, 0, &
                          1, npx-1, 1, npy-1, npz)
                                           call timing_off('COMM_TOTAL')
                  endif
               enddo

               !nested grid: receive from parent grids
               if (nested) then

                                           call timing_on('COMM_TOTAL')
                  call nested_grid_BC_save(parent_grid%w, nest_domain, &
                       ind_h, wt_h, 0, 0,  npx,  npy,  npz, 1, parent_grid%npx-1, 1, parent_grid%npy-1, &
                       var_east=w_east, &
                       var_west=w_west, &
                       var_north=w_north, &
                       var_south=w_south, ns=0, proc_in=.true.)
                                           call timing_off('COMM_TOTAL')
               endif
               
            else
                                           call timing_on('COMM_TOTAL')
               call nested_grid_BC_save(parent_grid%w, nest_domain, &
                 ind_h, wt_h, 0, 0, npx, npy, npz, isg, ieg, jsg, jeg, &
                 var_east=w_east, &
                 var_west=w_west, &
                 var_north=w_north, &
                 var_south=w_south, ns=nsponge, proc_in=ANY(pelist == gid)  )
                                           call timing_off('COMM_TOTAL')
            endif
            !delz not defined in halo

         end if
#endif

         !u
         if (concurrent) then

            !Coarse grid: send to child grids
            do p=1,size(child_grids)
               if (child_grids(p)) then
                                           call timing_on('COMM_TOTAL')
                  call nested_grid_BC_save(u, Atm(p)%nest_domain, 0, 1, &
                       1, npx-1, 1, npy-1, npz)
                                           call timing_off('COMM_TOTAL')
               endif
            enddo

            !nested grid: receive from parent grids
            if (nested) then

                                           call timing_on('COMM_TOTAL')
               call nested_grid_BC_save(parent_grid%u, nest_domain, &
                    ind_u, wt_u, 0, 1,  npx,  npy,  npz, 1, parent_grid%npx-1, 1, parent_grid%npy-1, &
                    var_east=u_east, &
                    var_west=u_west, &
                    var_north=u_north, &
                    var_south=u_south, ns=0, proc_in=.true.)
                                           call timing_off('COMM_TOTAL')
            endif

         else

                                           call timing_on('COMM_TOTAL')
            call nested_grid_BC_save(parent_grid%u, nest_domain, &
                 ind_u, wt_u, 0, 1, npx, npy, npz, isg, ieg, jsg, jeg, &
                 var_east=u_east, &
                 var_west=u_west, &
                 var_north=u_north, &
                 var_south=u_south, ns=nsponge, proc_in=ANY(pelist == gid)  )
                                           call timing_off('COMM_TOTAL')


         endif

         if (concurrent) then

            !Coarse grid: send to child grids
            do p=1,size(child_grids)
               if (child_grids(p)) then
                                           call timing_on('COMM_TOTAL')
                  call nested_grid_BC_save(vc, Atm(p)%nest_domain, 0, 1, &
                       1, npx-1, 1, npy-1, npz)
                                           call timing_off('COMM_TOTAL')
               endif
            enddo

            !nested grid: receive from parent grids
            if (nested) then
         
                                           call timing_on('COMM_TOTAL')
               call nested_grid_BC_save(parent_grid%vc, nest_domain, &
                    ind_u, wt_u, 0, 1,  npx,  npy,  npz, 1, parent_grid%npx-1, 1, parent_grid%npy-1, &
                    var_east=vc_east_t1, &
                    var_west=vc_west_t1, &
                    var_north=vc_north_t1, &
                    var_south=vc_south_t1, ns=0, proc_in=.true.)
                                           call timing_off('COMM_TOTAL')
            endif

         else
         
                                           call timing_on('COMM_TOTAL')
            call nested_grid_BC_save(parent_grid%vc, nest_domain, &
                 ind_u, wt_u, 0, 1, npx, npy, npz, isg, ieg, jsg, jeg, &
                 var_east=vc_east_t1, &
                 var_west=vc_west_t1, &
                 var_north=vc_north_t1, &
                 var_south=vc_south_t1, ns=nsponge, proc_in=ANY(pelist == gid)  )
                                           call timing_off('COMM_TOTAL')
         endif

         !v
         
         if (concurrent) then

            !Coarse grid: send to child grids
            do p=1,size(child_grids)
               if (child_grids(p)) then
                                           call timing_on('COMM_TOTAL')
                  call nested_grid_BC_save(v, Atm(p)%nest_domain, 1, 0, &
                       1, npx-1, 1, npy-1, npz)
                                           call timing_off('COMM_TOTAL')
               endif
            enddo

            !nested grid: receive from parent grids
            if (nested) then

                                           call timing_on('COMM_TOTAL')
               call nested_grid_BC_save(parent_grid%v, nest_domain, &
                    ind_v, wt_v, 1, 0,  npx,  npy,  npz, 1, parent_grid%npx-1, 1, parent_grid%npy-1, &
                    var_east=v_east, &
                    var_west=v_west, &
                    var_north=v_north, &
                    var_south=v_south, ns=0, proc_in=.true.)
                                           call timing_off('COMM_TOTAL')
            endif

         else
                                           call timing_on('COMM_TOTAL')
            call nested_grid_BC_save(parent_grid%v, nest_domain, &
              ind_v, wt_v, 1, 0, npx, npy, npz, isg, ieg, jsg, jeg, &
              var_east=v_east, &
              var_west=v_west, &
              var_north=v_north, &
              var_south=v_south, ns=nsponge, proc_in=ANY(pelist == gid)  )
                                           call timing_off('COMM_TOTAL')
         endif

         if (concurrent) then

            !Coarse grid: send to child grids
            do p=1,size(child_grids)
               if (child_grids(p)) then
                                           call timing_on('COMM_TOTAL')
                  call nested_grid_BC_save(uc, Atm(p)%nest_domain, 1, 0, &
                       1, npx-1, 1, npy-1, npz)
                                           call timing_off('COMM_TOTAL')
               endif
            enddo

            !nested grid: receive from parent grids
            if (nested) then

                                           call timing_on('COMM_TOTAL')
               call nested_grid_BC_save(parent_grid%uc, nest_domain, &
                    ind_v, wt_v, 1, 0,  npx,  npy,  npz, 1, parent_grid%npx-1, 1, parent_grid%npy-1, &
                    var_east=uc_east_t1, &
                    var_west=uc_west_t1, &
                    var_north=uc_north_t1, &
                    var_south=uc_south_t1, ns=0, proc_in=.true.)
                                           call timing_off('COMM_TOTAL')
            endif

         else

                                           call timing_on('COMM_TOTAL')
            call nested_grid_BC_save(parent_grid%uc, nest_domain, &
              ind_v, wt_v, 1, 0, npx, npy, npz, isg, ieg, jsg, jeg, &
              var_east=uc_east_t1, &
              var_west=uc_west_t1, &
              var_north=uc_north_t1, &
              var_south=uc_south_t1, ns=nsponge, proc_in=ANY(pelist == gid)  )
                                           call timing_off('COMM_TOTAL')
         endif


         if (first_step .and. concurrent) then
            if (nested) call set_BCs_t0(ncnst)
            first_step = .false.
         endif

         call mpp_sync_self

 end subroutine setup_nested_grid_BCs

 subroutine set_BCs_t0(ncnst)

   integer, intent(IN) :: ncnst

   if (concurrent) then

      h_east_t0 = h_east
      h_west_t0 = h_west
      h_north_t0 = h_north
      h_south_t0 = h_south
      if (ncnst > 0) then
         q_east_t0 = q_east
         q_west_t0 = q_west
         q_north_t0 = q_north
         q_south_t0 = q_south
      endif
#ifndef SW_DYNAMICS
      pt_east_t0 = pt_east
      pt_west_t0 = pt_west
      pt_north_t0 = pt_north
      pt_south_t0 = pt_south
      pt_east_t0 = pt_east
      pt_west_t0 = pt_west
      pt_north_t0 = pt_north
      pt_south_t0 = pt_south
      if (.not. hydrostatic) then
         w_east_t0 = w_east
         w_west_t0 = w_west
         w_north_t0 = w_north
         w_south_t0 = w_south
      endif
#endif
      u_east_t0 = u_east
      u_west_t0 = u_west
      u_north_t0 = u_north
      u_south_t0 = u_south
      v_east_t0 = v_east
      v_west_t0 = v_west
      v_north_t0 = v_north
      v_south_t0 = v_south

   endif

   vc_east_t0 = vc_east_t1
   vc_west_t0 = vc_west_t1
   vc_north_t0 = vc_north_t1
   vc_south_t0 = vc_south_t1
   uc_east_t0 = uc_east_t1
   uc_west_t0 = uc_west_t1
   uc_north_t0 = uc_north_t1
   uc_south_t0 = uc_south_t1


 end subroutine set_BCs_t0

 subroutine twoway_nest_update(npx, npy, npz, cp_air, zvir, ncnst, sphum,     &
                        u, v, w, omga, hydrostatic, pt, delp, q,   &
                        uc, vc, kappa, pkz, delz, ps, conv_theta_in)

    real, intent(IN) :: kappa, cp_air
    real, intent(IN) :: zvir

    integer, intent(IN) :: npx, npy, npz
    integer, intent(IN) :: ncnst, sphum
    logical, intent(IN) :: hydrostatic
    logical, intent(IN), OPTIONAL :: conv_theta_in

    real, intent(inout), dimension(isd:ied  ,jsd:jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(isd:ied+1,jsd:jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   isd:ied  ,jsd:jed  ,npz)  !  W (m/s)
    real, intent(inout) :: omga(isd:ied,jsd:jed,npz)      ! Vertical pressure velocity (pa/s)
    real, intent(inout) :: pt(  isd:ied  ,jsd:jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(isd:ied  ,jsd:jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   isd:ied  ,jsd:jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: uc(isd:ied+1,jsd:jed  ,npz) ! (uc,vc) mostly used as the C grid winds
    real, intent(inout) :: vc(isd:ied  ,jsd:jed+1,npz)

    real, intent(inout) :: pkz (is:ie,js:je,npz)             ! finite-volume mean pk
    real, intent(inout) :: delz(is:ie,js:je,npz)   ! delta-height (m); non-hydrostatic only
    real, intent(inout) :: ps  (isd:ied  ,jsd:jed)           ! Surface pressure (pascal)

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

    real :: qdp(   isd:ied  ,jsd:jed  ,npz)
    real, allocatable :: qdp_coarse(:,:,:)

    if (nestbctype > 1) upoff = 0

    !We update actual temperature, not theta.
    !If pt is actual temperature, set conv_theta to .false.
    if (present(conv_theta_in)) conv_theta = conv_theta_in

    rg = kappa*cp_air
    
    if (concurrent) then
       
       if ((.not. parent_proc) .and. (.not. child_proc)) return

    endif

    call mpp_get_data_domain( parent_grid%domain, &
         isd_p,  ied_p,  jsd_p,  jed_p  )
    call mpp_get_compute_domain( parent_grid%domain, &
         isc_p,  iec_p,  jsc_p,  jec_p  )

   r = refinement
   s = r/2 !rounds down (since r > 0)

   !delp/ps

   if (nestupdate /= 3 .and. nestupdate /= 4 &
        .and. nestupdate /= 5 .and. nestupdate /= 6) then


                                           call timing_on('COMM_TOTAL')
         call update_coarse_grid(parent_grid%delp, delp, nest_domain,&
              ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,:), &
              isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npz, 0, 0, &
              refinement, nestupdate, upoff, nsponge)
                                           call timing_off('COMM_TOTAL')

      call mpp_sync!self

#ifdef SW_DYNAMICS
      if (parent_proc) then
         do j=jsd_p,jed_p
            do i=isd_p,ied_p

               parent_grid%ps(i,j) = &
                    parent_grid%delp(i,j,1)/grav 

            end do
         end do
      endif
#endif

   end if

   if (nestupdate /= 3) then

      allocate(qdp_coarse(isd_p:ied_p,jsd_p:jed_p,npz))

      do n=1,ncnst

         if (child_proc) then
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

         if (parent_proc) then
            do k=1,npz
            do j=jsd_p,jed_p
            do i=isd_p,ied_p
               qdp_coarse(i,j,k) = parent_grid%q(i,j,k,n)*parent_grid%delp(i,j,k)
               !parent_grid%q(i,j,k,n) = qdp_coarse(i,j,k)/parent_grid%delp(i,j,k)
            enddo
            enddo
            enddo
         else
            qdp_coarse = 0.
         endif

                                           call timing_on('COMM_TOTAL')
!!$            call update_coarse_grid(parent_grid%q(:,:,:,n), &
!!$                 q(:,:,:,n), nest_domain, &
            call update_coarse_grid(qdp_coarse, qdp, nest_domain, &
                 ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,:), &
                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npz, 0, 0, &
                 refinement, nestupdate, upoff, nsponge)
                                           call timing_off('COMM_TOTAL')

               call mpp_sync!self

         if (parent_proc) then
            do k=1,npz
            do j=jsd_p,jed_p
            do i=isd_p,ied_p
               parent_grid%q(i,j,k,n) = qdp_coarse(i,j,k)/parent_grid%delp(i,j,k)
            enddo
            enddo
            enddo
         endif

      end do

      deallocate(qdp_coarse)

#ifndef SW_DYNAMICS
      if (conv_theta) then

         if (child_proc) then
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

                                           call timing_on('COMM_TOTAL')
            call update_coarse_grid(parent_grid%pt, &
                 t_nest, nest_domain, &
                 ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,:), &
                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npz, 0, 0, &
                 refinement, nestupdate, upoff, nsponge)
                                           call timing_off('COMM_TOTAL')
      else

                                           call timing_on('COMM_TOTAL')
            call update_coarse_grid(parent_grid%pt, &
              pt, nest_domain, &
              ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,:), &
              isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npz, 0, 0, &
              refinement, nestupdate, upoff, nsponge)
                                           call timing_off('COMM_TOTAL')

      endif !conv_theta

      call mpp_sync!self

      if (.not. hydrostatic) then

                                           call timing_on('COMM_TOTAL')
            call update_coarse_grid(parent_grid%w, w, nest_domain, &
                 ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,:), &
                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npz, 0, 0, &
                 refinement, nestupdate, upoff, nsponge)
            !Delz only defined in interior; will this work with the mpp nesting routines?
            call update_coarse_grid(parent_grid%delz, delz, nest_domain, &
                 ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,:), &
                 isd_p, ied_p, jsd_p, jed_p, isc, iec, jsc, jec, npz, 0, 0, &
                 refinement, nestupdate, upoff, nsponge)
                                           call timing_off('COMM_TOTAL')

         call mpp_sync!self

      end if
      

#endif
   end if !Nestupdate /= 3

                                           call timing_on('COMM_TOTAL')
      call update_coarse_grid(parent_grid%u, u, nest_domain, &
           ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,:), &
           isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npz, 0, 1, &
           refinement, nestupdate, upoff, nsponge)

      call update_coarse_grid(parent_grid%v, v, nest_domain, &
           ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,:), &
           isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npz, 1, 0, &
           refinement, nestupdate, upoff, nsponge)
                                           call timing_off('COMM_TOTAL')

   call mpp_sync!self


#ifndef SW_DYNAMICS
   if ((nestupdate == 5 .or. nestupdate == 6) .and. npz > 4) then

      !Use PS0 from nested grid, NOT the full delp. Also we assume the same number of levels on both grids.
      !PS0 should be initially set to be ps so that this routine does NOTHING outside of the update region

      !Re-compute nested (AND COARSE) grid ps

      allocate(ps0(isd_p:ied_p,jsd_p:jed_p))
      if (parent_proc) then

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

      if (child_proc) then

         ps = ptop
         do k=1,npz
            do j=jsd,jed
               do i=isd,ied
                  ps(i,j) = ps(i,j) + delp(i,j,k)
               end do
            end do
         end do
      endif

                                           call timing_on('COMM_TOTAL')
      call update_coarse_grid(ps0, ps, nest_domain,&
              ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,:), &
              isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, 0, 0, &
              refinement, nestupdate, upoff, nsponge)
                                           call timing_off('COMM_TOTAL')

      !!! The mpp version of update_coarse_grid does not return a consistent value of ps
      !!! across PEs, as it does not go into the haloes of a given coarse-grid PE. This
      !!! update_domains call takes care of the problem.

   if (parent_proc) then
                                                  call timing_on('COMM_TOTAL')
     call mpp_update_domains(ps0, parent_grid%domain, complete=.true.)
                                                 call timing_off('COMM_TOTAL')
   endif


      call mpp_sync!self

      if (parent_grid%tile == parent_tile) then 

         if (parent_proc) then

!!$         call update_remap_uv(npz, parent_grid%ak, parent_grid%bk, &
!!$              parent_grid%ps, &
!!$              parent_grid%u, &
!!$              parent_grid%v, npz, ps0, &
!!$              isc_p, iec_p, jsc_p, jec_p, isd_p, ied_p, jsd_p, jed_p)
         call update_remap_uv(npz, parent_grid%ak, parent_grid%bk, &
              parent_grid%ps(isd_p:ied_p,jsd_p:jed_p), &
              parent_grid%u(isd_p:ied_p,jsd_p:jed_p+1,1:npz), &
              parent_grid%v(isd_p:ied_p+1,jsd_p:jed_p,1:npz), npz, ps0(isd_p:ied_p,jsd_p:jed_p), &
              isc_p, iec_p, jsc_p, jec_p, isd_p, ied_p, jsd_p, jed_p)

         !comment out if statement to always remap theta instead of t in the remap-update.
         !(In LtE typically we use remap_t = .true.: remapping t is better (except in
         !idealized simulations with a background uniform theta) since near the top
         !boundary theta is exponential, which is hard to accurately interpolate with a spline
         if (.not. parent_grid%remap_t) then
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
              parent_grid%pt, parent_grid%q, npz, ps0, zvir, ncnst, &
              isc_p, iec_p, jsc_p, jec_p, isd_p, ied_p, jsd_p, jed_p)
         if (.not. parent_grid%remap_t) then
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

         endif !parent_proc

      end if 

      if (allocated(ps0)) deallocate(ps0)

   end if

#endif

   if (.not. concurrent) call setup_nested_BC_ucvc(current_Atm, npz, child_proc)

 end subroutine twoway_nest_update

!Do ua,va need to be converted back from lat-lon coords?
 subroutine before_twoway_nest_update(npx, npy, npz, ng, consv_te,               &
                        kappa, cp_air, zvir, ncnst,   &
                        u, v, w, delz, hydrostatic, pt, delp, q,   &
                        ps, pe, pk, peln, pkz, phis, ua, va, &
                        dry_mass, grid_number, mountain, make_nh)

    real, intent(IN) :: consv_te
    real, intent(IN) :: kappa, cp_air
    real, intent(IN) :: zvir

    integer, intent(IN) :: npx, npy, npz
    integer, intent(IN) :: ng
    integer, intent(IN) :: ncnst
    logical, intent(IN) :: hydrostatic

    real, intent(inout), dimension(isd:ied  ,jsd:jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(isd:ied+1,jsd:jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   isd:ied  ,jsd:jed  ,npz)  !  W (m/s)
    real, intent(inout) :: pt(  isd:ied  ,jsd:jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(isd:ied  ,jsd:jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   isd:ied  ,jsd:jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: delz(is:ie,js:je,npz)   ! delta-height (m); non-hydrostatic only

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout) :: ps  (isd:ied  ,jsd:jed)           ! Surface pressure (pascal)
    real, intent(inout) :: pe  (is-1:ie+1, npz+1,js-1:je+1)  ! edge pressure (pascal)
    real, intent(inout) :: pk  (is:ie,js:je, npz+1)          ! pe**cappa
    real, intent(inout) :: peln(is:ie,npz+1,js:je)           ! ln(pe)
    real, intent(inout) :: pkz (is:ie,js:je,npz)             ! finite-volume mean pk
    
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, intent(inout) :: phis(isd:ied,jsd:jed)       ! Surface geopotential (g*Z_surf)

    real, intent(inout), dimension(isd:ied ,jsd:jed ,npz):: ua, va

    real, intent(IN) :: dry_mass
    integer, intent(IN) :: grid_Number
    logical, intent(IN) :: mountain, make_nh

    real::   teq(is:ie,js:je)
    integer i, j, k, iq
    real rg

    integer :: sphum
    logical :: used

    sphum = get_tracer_index(MODEL_ATMOS, 'sphum')

   !Calculate cubed-sphere a-grid winds
   
    do k=1,npz
       call d2a_setup(u(isd,jsd,k), v(isd,jsd,k), ua(isd,jsd,k), va(isd,jsd,k), .true., &
            isd,ied,jsd,jed, is,ie,js,je, npx,npy, &
            grid_type, nested, cosa_s, rsin2)
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
        peln, phis, rsin2, cosa_s, zvir, cp_air,  rg, hlv, te_2d_coarse, &
        ua, va, teq, moist_phys, sphum, hydrostatic,idiag%id_te)

#endif

 end subroutine before_twoway_nest_update

 !NOTE: probably could also do the 'adjust dry mass' step here now
 subroutine after_twoway_nest_update(npz, nq,  ng, bdt, consv_te, fill,               &
                        reproduce_sum, kappa, cp_air, zvir, ks, ncnst,   &
                        u, v, w, delz, hydrostatic, pt, delp, q,   &
                        ps, pe, pk, peln, pkz, phis, omga, ua, va, uc, vc,          &
                        ak, bk, ze0, hybrid_z, dry_mass, adjust_dry_mass, grid_number, mountain, make_nh)

    real, intent(IN) :: bdt  ! Large time-step
    real, intent(IN) :: consv_te
    real, intent(IN) :: kappa, cp_air
    real, intent(IN) :: zvir

    integer, intent(IN) :: npz
    integer, intent(IN) :: nq             ! transported tracers
    integer, intent(IN) :: ng
    integer, intent(IN) :: ks
    integer, intent(IN) :: ncnst
    logical, intent(IN) :: fill
    logical, intent(IN) :: reproduce_sum
    logical, intent(IN) :: hydrostatic
    logical, intent(IN) :: hybrid_z       ! Using hybrid_z for remapping

    real, intent(inout), dimension(isd:ied  ,jsd:jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(isd:ied+1,jsd:jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   isd:ied  ,jsd:jed  ,npz)  !  W (m/s)
    real, intent(inout) :: pt(  isd:ied  ,jsd:jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(isd:ied  ,jsd:jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   isd:ied  ,jsd:jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: delz(is:ie,js:je,npz)   ! delta-height (m); non-hydrostatic only
    real, intent(inout) ::  ze0(is:ie,js:je,npz+1) ! height at edges (m); non-hydrostatic

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout) :: ps  (isd:ied  ,jsd:jed)           ! Surface pressure (pascal)
    real, intent(inout) :: pe  (is-1:ie+1, npz+1,js-1:je+1)  ! edge pressure (pascal)
    real, intent(inout) :: pk  (is:ie,js:je, npz+1)          ! pe**cappa
    real, intent(inout) :: peln(is:ie,npz+1,js:je)           ! ln(pe)
    real, intent(inout) :: pkz (is:ie,js:je,npz)             ! finite-volume mean pk
    
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, intent(inout) :: phis(isd:ied,jsd:jed)       ! Surface geopotential (g*Z_surf)
    real, intent(inout) :: omga(isd:ied,jsd:jed,npz)   ! Vertical pressure velocity (pa/s)
    real, intent(inout) :: uc(isd:ied+1,jsd:jed  ,npz) ! (uc,vc) mostly used as the C grid winds
    real, intent(inout) :: vc(isd:ied  ,jsd:jed+1,npz)

    real, intent(inout), dimension(isd:ied ,jsd:jed ,npz):: ua, va
    real, intent(in),    dimension(npz+1):: ak, bk

    real, intent(IN) :: dry_mass
    logical, intent(IN) :: adjust_dry_mass
    integer, intent(IN) :: grid_Number
    logical, intent(IN) :: mountain, make_nh

    real :: akap, tpe, rg
    integer:: kord_tracer(ncnst), cld_amt, iq
    real te_2d_coarse_after(is:ie,js:je)
    
    integer :: sphum, i, j, k
    real::   teq(is:ie,js:je)

    sphum = get_tracer_index(MODEL_ATMOS, 'sphum')

    cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
    akap  = kappa

#ifndef SW_DYNAMICS

   !To get coarse grid pkz, etc right after a two-way update so
   !that it is consistent across a restart:
   !(should only be called after doing such an update)

   call p_var(npz, is, ie, js, je, ptop, ptop_min,  &
        delp, delz, &
        pt, ps, &
        pe, peln,   &
        pk,   pkz, kappa, &
        q, ng, ncnst,  dry_mass,  &
        adjust_dry_mass,  mountain, &
        moist_phys,  hydrostatic, &
        nwat, make_nh)

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
        peln, phis, rsin2, cosa_s, zvir, cp_air,  rg, hlv, te_2d_coarse_after, &
        ua, va, teq, moist_phys, sphum, hydrostatic,idiag%id_te)

   te_2d_coarse = te_2d_coarse - te_2d_coarse_after
   tpe = g_sum(te_2d_coarse, is, ie, js, je, ng, area, 0)
   E_Flux_nest = tpe / (grav*bdt*4.*pi*radius**2)

#endif

 end subroutine after_twoway_nest_update

 !Routines for remapping (interpolated) nested-grid data to the coarse-grid's vertical coordinate.

 !This does not yet do anything for the tracers
 subroutine update_remap_tq( npz, ak,  bk,  ps, delp,  t,  q,  &
                      kmd, ps0, zvir, nq, &
                      is, ie, js, je, isd, ied, jsd, jed)
  integer, intent(in):: npz, kmd, nq
  real,    intent(in):: zvir
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in), dimension(isd:ied,jsd:jed):: ps0
  real,    intent(in), dimension(isd:ied,jsd:jed):: ps
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real,    intent(inout), dimension(isd:ied,jsd:jed,npz):: t
  real,    intent(inout), dimension(isd:ied,jsd:jed,npz,nq):: q
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
   call mappm(kmd, pn0, tp, npz, pn1, qn1, is,ie, 1, 8)

        do k=1,npz
           do i=is,ie
              t(i,j,k) = qn1(i,k)
           enddo
        enddo

5000 continue

 end subroutine update_remap_tq

 !remap_uv as-is remaps only a-grid velocities. A new routine has been written to handle staggered grids.
 subroutine update_remap_uv(npz, ak, bk, ps, u, v, kmd, ps0, &
                      is, ie, js, je, isd, ied, jsd, jed)
  integer, intent(in):: npz
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in):: ps(isd:ied,jsd:jed)
  real,    intent(inout), dimension(isd:ied,jsd:jed+1,npz):: u
  real,    intent(inout), dimension(isd:ied+1,jsd:jed,npz):: v
!
  integer, intent(in):: kmd
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
      call mappm(kmd, pe0(is:ie,:), qt(is:ie,:), npz, pe1(is:ie,:), qn1(is:ie,:), is,ie, -1, 8)
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
      call mappm(kmd, pe0(is:ie+1,:), qt(is:ie+1,:), npz, pe1(is:ie+1,:), qn1(is:ie+1,:), is,ie+1, -1, 8)
      do k=1,npz
         do i=is,ie+1
            v(i,j,k) = qn1(i,k)
         enddo
      enddo
   end do

 end subroutine update_remap_uv




end module fv_nesting_mod
