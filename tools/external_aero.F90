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

! =======================================================================
! this module is designed to read 12 months climatology aerosol and
! interpolate to daily aerosol
! developer: linjiong zhou
! =======================================================================

module external_aero_mod

	use fms_mod, only: mpp_error, FATAL
        use fms2_io_mod, only: file_exists
	use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_npes, mpp_get_current_pelist
	use time_manager_mod, only: time_type
	use fv_mapz_mod, only: map1_q2
	use fv_fill_mod, only: fillz

	public :: load_aero, read_aero, clean_aero

	! MERRA2 aerosol: # month = 12, # vertical layer = 72
	integer :: nmon = 12, nlev = 72
	integer :: id_aero, id_aero_now

	! share arrays for time and level interpolation
	real, allocatable, dimension(:,:,:) :: aero_ps
	real, allocatable, dimension(:,:,:,:) :: aero_p
	real, allocatable, dimension(:,:,:,:) :: aero_pe
	real, allocatable, dimension(:,:,:,:) :: aero_dp
	real, allocatable, dimension(:,:,:,:) :: aerosol

contains

! =======================================================================
! load aerosol 12 months climatological dataset

subroutine load_aero(Atm, Time)

        use fms2_io_mod, only: FmsNetcdfDomainFile_t, open_file, close_file, &
                               register_restart_field, register_axis, &
                               read_restart, get_variable_dimension_names, &
                               get_dimension_size, close_file
	use fv_arrays_mod, only: fv_atmos_type
	use diag_manager_mod, only: register_static_field, register_diag_field

	implicit none

	type(time_type), intent(in) :: Time
	type(fv_atmos_type), intent(in), target :: Atm
        type(FmsNetcdfDomainFile_t) :: aero_restart

	integer :: k
	integer :: is, ie, js, je
        character(len=8), dimension(4) :: dim_names_4d
        character(len=8), dimension(3) :: dim_names_3d
        integer, dimension(2) :: dim_size

	real, allocatable, dimension(:,:,:,:) :: aero_lndp

	character(len=64) :: file_name = "MERRA2_400.inst3_3d_aer_Nv.climatology.nc"

	is = Atm%bd%is
	ie = Atm%bd%ie
	js = Atm%bd%js
	je = Atm%bd%je

	if (mpp_pe() .eq. mpp_root_pe()) then
		write(*,*) "aerosol 12 months climatological dataset is used for forecast."
	endif

        ! -----------------------------------------------------------------------
        ! load aerosol data

        if (open_file(aero_restart, 'INPUT/'//trim(file_name), "read", Atm%domain, &
                      & is_restart=.true., dont_add_res_to_filename=.true.)) then

		! allocate share arrays
		if (.not. allocated(aero_ps)) allocate(aero_ps(is:ie,js:je,nmon))
		if (.not. allocated(aero_p)) allocate(aero_p(is:ie,js:je,nlev,nmon))
		if (.not. allocated(aero_pe)) allocate(aero_pe(is:ie,js:je,nlev+1,nmon))
		if (.not. allocated(aero_dp)) allocate(aero_dp(is:ie,js:je,nlev,nmon))
		if (.not. allocated(aerosol)) allocate(aerosol(is:ie,js:je,nlev,nmon))

                ! read in restart files
                call get_variable_dimension_names(aero_restart, "PS", dim_names_3d)
                call get_variable_dimension_names(aero_restart, "DELP", dim_names_4d)
                call get_dimension_size(aero_restart, dim_names_4d(3), dim_size(1))
                call get_dimension_size(aero_restart, dim_names_4d(4), dim_size(2))
                call register_axis(aero_restart, dim_names_4d(1), "x")
                call register_axis(aero_restart, dim_names_4d(2), "y")
                call register_axis(aero_restart, dim_names_4d(3), dim_size(1))
                call register_axis(aero_restart, dim_names_4d(4), dim_size(2))
                call register_restart_field(aero_restart,"PS",&
                 aero_ps, dim_names_3d)
                call register_restart_field(aero_restart,"DELP",&
                aero_dp, dim_names_4d)
                call register_restart_field(aero_restart,"SO4",&
                aerosol, dim_names_4d)
                call read_restart(aero_restart)
                call close_file(aero_restart)
        else

		! stop when aerosol does not exist
		call mpp_error("external_aero_mod",&
			"file: "//trim(file_name)//" does not exist.",FATAL)

	endif

	! -----------------------------------------------------------------------
	! calculate layer mean pressure

	! allocate local array
	if (.not. allocated(aero_lndp)) allocate(aero_lndp(is:ie,js:je,nlev,nmon))

	! calcuate edge pressure
	aero_p = -999.9
	aero_pe(:,:,nlev+1,:) = aero_ps
	do k = nlev, 1, -1
		aero_pe(:,:,k,:) = aero_pe(:,:,k+1,:) - aero_dp(:,:,k,:)
	enddo

	! stop when minimum value is less and equal to zero
	if (minval(aero_pe) .le. 0.0) then
		call mpp_error("external_aero_mod","aero_pe has value <= 0.",FATAL)
	endif

	! calcuate layer mean pressure
	do k = 1, nlev
		aero_lndp(:,:,k,:) = log(aero_pe(:,:,k+1,:)) - log(aero_pe(:,:,k,:))
	enddo
	aero_p = aero_dp / aero_lndp

	! stop when minimum value is less and equal to zero
	if (minval(aero_p) .le. 0.0) then
		call mpp_error("external_aero_mod","aero_p has value <= 0.",FATAL)
	endif

	! deallocate local array
	if (allocated(aero_lndp)) deallocate(aero_lndp)

	! -----------------------------------------------------------------------
	! register for diagnostic output

	id_aero = register_static_field('dynamics','aero_ann',&
		Atm%atmos_axes(1:2),'none','none')
	id_aero_now= register_diag_field('dynamics','aero_now',&
		Atm%atmos_axes(1:2),Time,'none','none')

end subroutine load_aero

! =======================================================================
! read aerosol climatological dataset

subroutine read_aero(is, ie, js, je, npz, nq, Time, pe, peln, qa, kord_tr, fill)

	use constants_mod, only: grav
	use diag_manager_mod, only: send_data
	use time_manager_mod, only: get_date, set_date, get_time, operator(-)
	use tracer_manager_mod, only: get_tracer_index
	use field_manager_mod, only: MODEL_ATMOS

	implicit none

	type(time_type), intent(in) :: Time
	type(time_type) :: Time_before
	type(time_type) :: Time_after

	integer :: i, j, k, n
	integer, intent(in) :: is, ie, js, je, npz, nq, kord_tr
	integer :: year, month, day, hour, minute, second
	integer :: seconds, days01, days21, month1, month2
	integer :: aero_id

	real, dimension(is:ie,js:je,npz,nq), intent(inout) :: qa
	real, dimension(is:ie,npz+1,js:je), intent(in) :: pe, peln

	real, allocatable, dimension(:,:) :: vi_aero
	real, allocatable, dimension(:,:) :: vi_aero_now
	real, allocatable, dimension(:,:,:) :: aero_now_a
	real, allocatable, dimension(:,:,:) :: aero_now_p
	real, allocatable, dimension(:,:,:) :: aero_now_pe
	real, allocatable, dimension(:,:,:) :: aero_now_dp
	real, allocatable, dimension(:,:,:) :: pm

	logical :: used, use_fv3_interp = .true.
	logical, intent (in) :: fill

	! -----------------------------------------------------------------------
	! diagnostic output of annual mean vertical integral aerosol

	if (id_aero > 0) then

		! allocate local array
		if (.not. allocated(vi_aero)) allocate(vi_aero(is:ie,js:je))
  
		! calcualte annual mean vertical intergral aerosol
		vi_aero = 0.0
		do n = 1, nmon
			do k = 1, nlev
				vi_aero = vi_aero + aerosol(:,:,k,n) * aero_dp(:,:,k,n)
			enddo
		enddo
		vi_aero = vi_aero / nmon / grav * 1.e6
  
		! diagnostic output
		used = send_data(id_aero,vi_aero,Time)
  
		! deallocate local array
		if (allocated(vi_aero)) deallocate(vi_aero)

	endif

	! -----------------------------------------------------------------------
	! linearly interpolate monthly aerosol to today

	! allocate local array
	if (.not. allocated(aero_now_a)) allocate(aero_now_a(is:ie,js:je,nlev))
	if (.not. allocated(aero_now_p)) allocate(aero_now_p(is:ie,js:je,nlev))
	if (.not. allocated(aero_now_pe)) allocate(aero_now_pe(is:ie,js:je,nlev+1))

	! get current date information
	call get_date(Time, year, month, day, hour, minute, second)

	! get previous day 15 and next day 15 time
	if (day .ge. 15) then
		Time_before = set_date(year, month, 15, 0, 0, 0)
		if (month .eq. 12) then
			Time_after = set_date(year+1, 1, 15, 0, 0, 0)
		else
			Time_after = set_date(year, month+1, 15, 0, 0, 0)
		endif
	else
		if (month .eq. 1) then
			Time_before = set_date(year-1, 12, 15, 0, 0, 0)
		else
			Time_before = set_date(year, month-1, 15, 0, 0, 0)
		endif
		Time_after = set_date(year, month, 15, 0, 0, 0)
	endif

	! get day difference between current day and previous day 15,
	! and between next day 15 and previous day 15
	call get_time(Time - Time_before, seconds, days01)
	call get_time(Time_after - Time_before, seconds, days21)
	call get_date(Time_before, year, month1, day, hour, minute, second)
	call get_date(Time_after, year, month2, day, hour, minute, second)

	! get aerosol for current date
	aero_now_a = aerosol(:,:,:,month2) - aerosol(:,:,:,month1)
	aero_now_a = 1.0 * days01 / days21 * aero_now_a + aerosol(:,:,:,month1)
	aero_now_p = aero_p(:,:,:,month2) - aero_p(:,:,:,month1)
	aero_now_p = 1.0 * days01 / days21 * aero_now_p + aero_p(:,:,:,month1)
	aero_now_pe = aero_pe(:,:,:,month2) - aero_pe(:,:,:,month1)
	aero_now_pe = 1.0 * days01 / days21 * aero_now_pe + aero_pe(:,:,:,month1)

	! -----------------------------------------------------------------------
	! diagnostic output of current vertical integral aerosol

	if (id_aero_now > 0) then

		! allocate local array
		if (.not. allocated(vi_aero_now)) allocate(vi_aero_now(is:ie,js:je))
		if (.not. allocated(aero_now_dp)) allocate(aero_now_dp(is:ie,js:je,nlev))
  
		! get pressure thickness for current date
		aero_now_dp = aero_dp(:,:,:,month2) - aero_dp(:,:,:,month1)
		aero_now_dp = 1.0 * days01 / days21 * aero_now_dp + aero_dp(:,:,:,month1)

		! calcualte annual mean vertical intergral aerosol
		vi_aero_now = 0.0
		do k = 1, nlev
			vi_aero_now = vi_aero_now + aero_now_a(:,:,k) * aero_now_dp(:,:,k)
		enddo
		vi_aero_now = vi_aero_now / grav * 1.e6
  
		! diagnostic output
		used = send_data(id_aero_now,vi_aero_now,Time)
  
		! deallocate local array
		if (allocated(vi_aero_now)) deallocate(vi_aero_now)
		if (allocated(aero_now_dp)) deallocate(aero_now_dp)

	endif

	! -----------------------------------------------------------------------
	! vertically interpolate aeorosol

	! allocate local array
	if (.not. allocated(pm)) allocate(pm(is:ie,js:je,npz))

	! calculate layer mean pressure
	do k = 1, npz
		pm(:,:,k) = (pe(:,k+1,:) - pe(:,k,:)) / (peln(:,k+1,:) - peln(:,k,:))
	enddo

	! stop when minimum value is less and equal to zero
	if (minval(pm) .le. 0.0) then
		call mpp_error("external_aero_mod","pm has value <= 0.",FATAL)
	endif

	! get aerosol tracer id
	aero_id = get_tracer_index(MODEL_ATMOS, 'aerosol')

	! vertically interpolation
	if (use_fv3_interp) then
		do j = js, je
			call map1_q2 (nlev, aero_now_pe (is:ie, j, :), aero_now_a (is:ie, js:je, :), &
				npz, pe (is:ie, :, j), qa (is:ie, j, :, aero_id), &
				pe (is:ie, 2:npz+1, j) - pe (is:ie, 1:npz, j), &
				is, ie, 0, kord_tr, j, is, ie, js, je, 0.)
			if (fill) call fillz (ie-is+1, npz, 1, qa (is:ie, j, :, aero_id), &
				pe (is:ie, 2:npz+1, j) - pe (is:ie, 1:npz, j))
		enddo
	else
		do j = js, je
			do i = is, ie
				do k = 1, npz
					if (pm(i,j,k) .lt. aero_now_p(i,j,1)) then
						qa(i,j,k,aero_id) = aero_now_a(i,j,1)
						!qa(i,j,k,aero_id) = aero_now_a(i,j,1) + &
						!	(log(pm(i,j,k)) - log(aero_now_p(i,j,1))) / &
						!	(log(aero_now_p(i,j,2)) - log(aero_now_p(i,j,1))) * &
						!	(aero_now_a(i,j,2) - aero_now_a(i,j,1))
					else if (pm(i,j,k) .ge. aero_now_p(i,j,nlev)) then
						qa(i,j,k,aero_id) = aero_now_a(i,j,nlev)
						!qa(i,j,k,aero_id) = aero_now_a(i,j,nlev-1) + &
						!	(log(pm(i,j,k)) - log(aero_now_p(i,j,nlev-1))) / &
						!	(log(aero_now_p(i,j,nlev)) - log(aero_now_p(i,j,nlev-1))) * &
						!	(aero_now_a(i,j,nlev) - aero_now_a(i,j,nlev-1))
					else
						do n = 1, nlev-1
							if (pm(i,j,k) .ge. aero_now_p(i,j,n) .and. &
								pm(i,j,k) .lt. aero_now_p(i,j,n+1)) then
								qa(i,j,k,aero_id) = aero_now_a(i,j,n) + &
									(log(pm(i,j,k)) - log(aero_now_p(i,j,n))) / &
									(log(aero_now_p(i,j,n+1)) - log(aero_now_p(i,j,n))) * &
									(aero_now_a(i,j,n+1) - aero_now_a(i,j,n))
							endif
						enddo
					endif
				enddo
			enddo
		enddo
	endif

	! deallocate local array
	if (allocated(pm)) deallocate(pm)

	! -----------------------------------------------------------------------
	! deallocate local array

	if (allocated(aero_now_a)) deallocate(aero_now_a)
	if (allocated(aero_now_p)) deallocate(aero_now_p)

end subroutine read_aero

! =======================================================================
! clean aerosol climatological dataset

subroutine clean_aero()

	implicit none

	if (allocated(aero_ps)) deallocate(aero_ps)
	if (allocated(aero_p)) deallocate(aero_p)
	if (allocated(aero_pe)) deallocate(aero_pe)
	if (allocated(aero_dp)) deallocate(aero_dp)
	if (allocated(aerosol)) deallocate(aerosol)

end subroutine clean_aero

end module external_aero_mod
