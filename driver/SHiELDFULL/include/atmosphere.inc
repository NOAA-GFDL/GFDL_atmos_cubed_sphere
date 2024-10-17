 subroutine ATMOSPHERE_GRID_BDRY_ (blon, blat, global)
!---------------------------------------------------------------
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid
!---------------------------------------------------------------
    real(ATMOSPHERE_KIND_),    intent(out) :: blon(:,:), blat(:,:)   ! Unit: radian
    logical, intent(in), optional :: global
! Local data:
    integer i,j

    if( PRESENT(global) ) then
      if (global) call mpp_error(FATAL, '==> global grid is no longer available &
                               & in the Cubed Sphere')
    endif

    if (ATMOSPHERE_KIND_ .eq. r8_kind) then
      do j=jsc,jec+1
        do i=isc,iec+1
          blon(i-isc+1,j-jsc+1) = Atm(mygrid)%gridstruct%grid_64(i,j,1)
          blat(i-isc+1,j-jsc+1) = Atm(mygrid)%gridstruct%grid_64(i,j,2)
        enddo
      end do
    else
      do j=jsc,jec+1
        do i=isc,iec+1
          blon(i-isc+1,j-jsc+1) = Atm(mygrid)%gridstruct%grid(i,j,1)
          blat(i-isc+1,j-jsc+1) = Atm(mygrid)%gridstruct%grid(i,j,2)
        enddo
      end do
    endif
 end subroutine ATMOSPHERE_GRID_BDRY_

 subroutine ATMOSPHERE_PREF_ (p_ref)
   real(ATMOSPHERE_KIND_), dimension(:,:), intent(inout) :: p_ref

   p_ref = _DBL_(_RL_(pref))

 end subroutine ATMOSPHERE_PREF_

  subroutine ATMOSPHERE_CELL_AREA_  (area_out)
   real(ATMOSPHERE_KIND_), dimension(:,:),  intent(out)          :: area_out

   area_out(1:iec-isc+1, 1:jec-jsc+1) =  _DBL_(_RL_(Atm(mygrid)%gridstruct%area (isc:iec,jsc:jec)))

 end subroutine ATMOSPHERE_CELL_AREA_

 subroutine GET_BOTTOM_MASS_ ( t_bot, tr_bot, p_bot, z_bot, p_surf, slp )
  !--------------------------------------------------------------
  ! returns temp, sphum, pres, height at the lowest model level
  ! and surface pressure
  !--------------------------------------------------------------
     real(ATMOSPHERE_KIND_), intent(out), dimension(isc:iec,jsc:jec):: t_bot, p_bot, z_bot, p_surf
     real(ATMOSPHERE_KIND_), intent(out), optional, dimension(isc:iec,jsc:jec):: slp
     real(ATMOSPHERE_KIND_), intent(out), dimension(isc:iec,jsc:jec,nq):: tr_bot
     integer :: i, j, m, k, kr
     real(ATMOSPHERE_KIND_)    :: rrg, sigtop, sigbot
     real(ATMOSPHERE_KIND_), dimension(isc:iec,jsc:jec) :: tref
     real(ATMOSPHERE_KIND_), parameter :: tlaps = 6.5e-3

     rrg  = _DBL_(_RL_(rdgas / grav))

     do j=jsc,jec
        do i=isc,iec
           p_surf(i,j) = _DBL_(_RL_(Atm(mygrid)%ps(i,j)))
           t_bot(i,j) = _DBL_(_RL_(Atm(mygrid)%pt(i,j,npz)))
           p_bot(i,j) = _DBL_(_RL_(Atm(mygrid)%delp(i,j,npz)/(Atm(mygrid)%peln(i,npz+1,j)-Atm(mygrid)%peln(i,npz,j))))
           z_bot(i,j) = rrg*t_bot(i,j)*_DBL_(_RL_((1.+zvir*Atm(mygrid)%q(i,j,npz,sphum)))) *  &
                        _DBL_(_RL_((1. - Atm(mygrid)%pe(i,npz,j)/p_bot(i,j))))
        enddo
     enddo

     if ( present(slp) ) then
       ! determine 0.8 sigma reference level
       sigtop = _DBL_(_RL_(Atm(mygrid)%ak(1)/pstd_mks+Atm(mygrid)%bk(1)))
       do k = 1, npz
          sigbot = _DBL_(_RL_(Atm(mygrid)%ak(k+1)/pstd_mks+Atm(mygrid)%bk(k+1)))
          if (sigbot+sigtop > 1.6) then
             kr = k
             exit
          endif
          sigtop = sigbot
       enddo
       do j=jsc,jec
          do i=isc,iec
             ! sea level pressure
             tref(i,j) = _DBL_(_RL_(Atm(mygrid)%pt(i,j,kr) * (Atm(mygrid)%delp(i,j,kr)/ &
                              ((Atm(mygrid)%peln(i,kr+1,j)-Atm(mygrid)%peln(i,kr,j))*Atm(mygrid)%ps(i,j)))**(-rrg*tlaps)))
             slp(i,j) = _DBL_(_RL_(Atm(mygrid)%ps(i,j)*(1.+tlaps*Atm(mygrid)%phis(i,j)/(real(tref(i,j))*grav))**(1./(rrg*tlaps))))
          enddo
       enddo
     endif

  ! Copy tracers
     do m=1,nq
        do j=jsc,jec
           do i=isc,iec
              tr_bot(i,j,m) = _DBL_(_RL_(Atm(mygrid)%q(i,j,npz,m)))
           enddo
        enddo
     enddo

   end subroutine GET_BOTTOM_MASS_

   subroutine GET_BOTTOM_WIND_ ( u_bot, v_bot )
!-----------------------------------------------------------
! returns u and v on the mass grid at the lowest model level
!-----------------------------------------------------------
   real(ATMOSPHERE_KIND_), intent(out), dimension(isc:iec,jsc:jec):: u_bot, v_bot
   integer i, j

   do j=jsc,jec
      do i=isc,iec
         u_bot(i,j) = _DBL_(_RL_(Atm(mygrid)%u_srf(i,j)))
         v_bot(i,j) = _DBL_(_RL_(Atm(mygrid)%v_srf(i,j)))
      enddo
   enddo

 end subroutine GET_BOTTOM_WIND_

 subroutine GET_STOCK_PE_(index, value)
   integer, intent(in) :: index
   real(ATMOSPHERE_KIND_),   intent(out) :: value

#ifdef USE_STOCK
   include 'stock.inc'
#endif

   real(ATMOSPHERE_KIND_) wm(isc:iec,jsc:jec)
   integer i,j,k

   select case (index)

#ifdef USE_STOCK
   case (ISTOCK_WATER)
#else
   case (1)
#endif

!----------------------
! Perform vertical sum:
!----------------------
     wm = 0.
     do j=jsc,jec
        do k=1,npz
           do i=isc,iec
! Warning: the following works only with AM2 physics: water vapor; cloud water, cloud ice.
              wm(i,j) = wm(i,j) + _DBL_(_RL_(Atm(mygrid)%delp(i,j,k) * ( Atm(mygrid)%q(i,j,k,sphum)   +  &
                                                              Atm(mygrid)%q(i,j,k,liq_wat) +  &
                                                              Atm(mygrid)%q(i,j,k,ice_wat) )))
           enddo
        enddo
     enddo

!----------------------
! Horizontal sum:
!----------------------
     value = 0.
     do j=jsc,jec
        do i=isc,iec
           value = value + wm(i,j)*_DBL_(_RL_(Atm(mygrid)%gridstruct%area(i,j)))
        enddo
     enddo
     value = value/_DBL_(_RL_(grav))

   case default
     value = 0.0
   end select

 end subroutine GET_STOCK_PE_
