
module w_forcing_mod

  use fv_arrays_mod,   only: fv_grid_type, fv_flags_type, fv_grid_bounds_type, R_GRID
  use mpp_domains_mod, only: mpp_update_domains, domain2d
  use mpp_mod,         only: mpp_error, FATAL, mpp_root_pe, mpp_broadcast, mpp_sum, mpp_sync
  use fv_mp_mod,       only: is_master
  implicit none
  public

  !settings
  integer :: w_forcing_type = 101
  real :: w_forcing_L = 40000. !m
  real :: w_forcing_R = 2000.  !m --- ?!?
  real :: w_forcing_D = 4000.  !m, depth
  real :: w_forcing_H = 0.
  real :: w_forcing_start = 0.0   !s
  real :: w_forcing_end  = -1. !2400.0 !s
  real :: w_forcing_a = 2.0 !acceleration, m/s**2

  real :: w_forcing_Divg = 3.75e-6 !1/s
  real :: w_forcing_tau = 3600. !s

  !saved data
  real :: w_forcing_i0
  real :: w_forcing_j0


contains

 subroutine init_w_forcing(bd, npx, npy, npz, grid_type, agrid, flagstruct)!, wft)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real  , intent(IN)      :: agrid(bd%isd:bd%ied, bd%jsd:bd%jed)
   integer,intent(IN)      :: npx, npy, npz, grid_type!, wft
   type(fv_flags_type), target, intent(IN) :: flagstruct

   !w_forcing_type = wft

   if (grid_type == 4) then

      select case (w_forcing_type)
      case(1)  ! half-ellipse acceleration (Ziegler et al., 2010; Prein et al. 2021)
         w_forcing_i0 = real(npx-1)*0.5
         w_forcing_j0 = real(npy-1)*0.5
      case default

      end select

   endif

   if (is_master()) print*, ' CALLING INIT_W_FORCING ', w_forcing_type, w_forcing_i0, w_forcing_j0

 end subroutine init_w_forcing

 subroutine do_w_forcing(bd, npx, npy, npz, w, delz, phis, grid_type, agrid, domain, flagstruct, dt, time)

   implicit none

   type(fv_grid_bounds_type), intent(IN) :: bd
   real  , intent(INOUT)   ::     w(bd%isd:, bd%jsd:,1:)
   real  , intent(IN)      ::  delz(bd%is:  , bd%js:  ,1:)
   real  , intent(IN)      ::  phis(bd%isd:bd%ied, bd%jsd:bd%jed)
   real  , intent(IN)      :: agrid(bd%isd:bd%ied, bd%jsd:bd%jed,2)
   integer,intent(IN)      :: npx, npy, npz, grid_type
   real  , intent(IN)      :: dt, time
   type(fv_flags_type), target, intent(IN) :: flagstruct
   type(domain2d), intent(INOUT) :: domain

   real :: Htop(bd%is:bd%ie,bd%js:bd%je) !height at the top of the current layer
   real :: rad,radm1,ht,xL,wls,forc,dttau,lev

   integer :: i,j,k

   integer :: is, ie, js, je
   integer :: isd, ied, jsd, jed

   logical, SAVE :: first_time = .true.

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed

   if (grid_type < 4) then
      call mpp_error(FATAL, "Not implemented for grid_type < 4 yet.")
   endif

   if (w_forcing_end > 0) then
      if (time < w_forcing_start .or. time > w_forcing_end) return
   endif

   if (first_time .and. is_master()) print*, ' CALLING DO_W_FORCING '

   if (grid_type == 4) then

      select case (w_forcing_type)
      case(1)


         do j=js,je
            do i=is,ie
               Htop(i,j) = phis(i,j) !Height above MSL
            enddo
         enddo
         do k=npz,1,-1
            do j=js,je
               do i=is,ie
                  Htop(i,j) = Htop(i,j) - delz(i,j,k)

                  xL = abs(i-w_forcing_i0)*flagstruct%dx_const
                  if (xL <= w_forcing_L) then
                     rad = (j-w_forcing_j0)*flagstruct%dx_const
                     rad = rad*rad/(w_forcing_R*w_forcing_R)
                     ht = Htop(i,j) + 0.5*delz(i,j,k) - w_forcing_H
                     rad = rad + ht*ht/(w_forcing_D*w_forcing_D)
                     radm1 = max(1.-sqrt(rad),0.)
                     w(i,j,k) = w(i,j,k) + w_forcing_a*radm1*radm1*dt
                  endif

               enddo
            enddo
         enddo

      case(101)
         !PBL simulations with specified divergence
         !Nudging domain to w = Dz
         !do not apply in sponge layer

         dttau=dt/w_forcing_tau
         forc = 1./(1.+dttau)
         do j=js,je
            do i=is,ie
               Htop(i,j) = -delz(i,j,npz)*0.5
               wls = -w_forcing_Divg*Htop(i,j)
               w(i,j,npz) = (w(i,j,npz) + dttau*wls)*forc
            enddo
         enddo
         do k=npz-1,3,-1
            do j=js,je
               do i=is,ie
                  Htop(i,j) = Htop(i,j) - 0.5*(delz(i,j,k-1)+delz(i,j,k))
                  wls = -w_forcing_Divg*Htop(i,j)
                  w(i,j,k) = (w(i,j,k) + dttau*wls)*forc
               enddo
            enddo
         enddo

         if (first_time .and. is_master()) then
            i=is
            j=js
            lev=-delz(i,j,npz)*0.5
            wls = -w_forcing_Divg*lev
            print*, npz, wls, w(i,j,npz), dttau
            do k=npz,3,-1
               lev = lev - 0.5*(delz(i,j,k-1)+delz(i,j,k))
               wls = -w_forcing_divg*lev
               print*, k, wls, w(i,j,k)
            enddo
         endif

      case default
         call mpp_error(FATAL, " Value of w_forcing_type not implemented.")

      end select

   end if

   call mpp_update_domains(w, domain)

   first_time = .false.

 end subroutine do_w_forcing


end module w_forcing_mod
