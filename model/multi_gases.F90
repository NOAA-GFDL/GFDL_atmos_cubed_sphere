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

!>@brief The module 'multi_gases' peforms multi constitutents computations.
!>@author H.-M. H. Juang, NOAA/NWS/NCEP/EMC

module multi_gases_mod

! Modules Included:
! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>rdgas, cp_air</td>
!   </tr>
! </table>

      use constants_mod,     only: rdgas, rvgas, cp_air
      use     fv_mp_mod,     only: is_master
      use mpp_mod,           only: stdlog, input_nml_file
      use fms_mod,           only: check_nml_error


      implicit none
      integer num_gas
      integer ind_gas
      integer num_wat
      integer sphum, sphump1
      real, allocatable :: ri(:)
      real, allocatable :: cpi(:)
      real, allocatable :: vir(:)
      real, allocatable :: vicp(:)
      real, allocatable :: vicv(:)

      private num_wat, sphum, sphump1
      public vir, vicp, vicv, ind_gas, num_gas
      public multi_gases_init, read_namelist_multi_gases_nml
      public virq
      public virq_max
      public virqd
      public vicpqd
      public virq_nodq
      public virq_qpz
      public vicpqd_qpz
      public vicvqd
      public vicvqd_qpz

      CONTAINS
! --------------------------------------------------------
      subroutine multi_gases_init(ngas, nwat)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: vir(i): ri/rdgas - r0/rdgas
!        vir(0): r0/rdgas
!        vicp(i): cpi/cp_air - cp0/cp_air
!        vicp(0): cp0/cp_air
!        cv_air = cp_air - rdgas
!        vicv(i): cvi/cv_air - cv0/cv_air
!        vicv(0): cv0/cv_air
!--------------------------------------------
      integer, intent(in):: ngas, nwat
! Local:
      integer n
      real   cvi(0:ngas)
      real   cv_air
      logical :: default_gas=.false.

      sphum = 1
      sphump1 = sphum+1
      num_wat = nwat
      ind_gas = num_wat+1
      do n=0,ngas
        if( ri(n).ne.0.0 .or. cpi(n).ne.0.0 ) num_gas=n
      enddo
      if ( num_gas.eq.1 ) default_gas=.true.
      allocate( vir (0:num_gas) )
      allocate( vicp(0:num_gas) )
      allocate( vicv(0:num_gas) )

      cv_air = cp_air - rdgas
      do n=0,num_gas
        cvi(n) = cpi(n) - ri(n)
      enddo

      vir (0) =  ri(0)/rdgas
      vicp(0) = cpi(0)/cp_air
      vicv(0) = cvi(0)/cv_air
      if( default_gas ) then
        vir (0) = 1.0
        vicp(0) = 1.0
        vicv(0) = 1.0
      endif
      do n=1,num_gas
        vir(n) = 0.0
        if(  ri(n).gt.0.0 ) vir (n) =  ri(n)/rdgas -  vir (0)
        vicp(n) = 0.0
        if( cpi(n).gt.0.0 ) vicp(n) = cpi(n)/cp_air - vicp(0)
        vicv(n) = 0.0
        if( cvi(n).gt.0.0 ) vicv(n) = cvi(n)/cv_air - vicv(0)
      enddo

      if( is_master() ) then
        write(*,*) ' multi_gases_init with ind_gas=',ind_gas
        write(*,*) ' multi_gases_init with num_gas=',num_gas
        write(*,*) ' multi_gases_init with vir =',vir
        write(*,*) ' multi_gases_init with vicp=',vicp
        write(*,*) ' multi_gases_init with vicv=',vicv
      endif

      return
      end subroutine multi_gases_init
      subroutine read_namelist_multi_gases_nml(nml_filename,ncnst,nwat)

        character(*), intent(IN) :: nml_filename
        integer, intent(IN) :: ncnst, nwat
        integer :: ierr, f_unit, unit, ios

        namelist /multi_gases_nml/ ri,cpi

       unit = stdlog()

       allocate (ri(0:ncnst))
       allocate (cpi(0:ncnst))

       ri     =    0.0
       cpi    =    0.0
       ri(0)  = rdgas
       ri(1)  = rvgas
       cpi(0) = cp_air
       cpi(1) = 4*cp_air
#ifdef INTERNAL_FILE_NML

      ! Read multi_gases namelist
        read (input_nml_file,multi_gases_nml,iostat=ios)
        ierr = check_nml_error(ios,'multi_gases_nml')
                                                             
#else
      ! Read multi_gases namelist
        f_unit = open_namelist_file(nml_filename)

        rewind (f_unit)
        read (f_unit,multi_gases_nml,iostat=ios)
        ierr = check_nml_error(ios,'multi_gases_nml')
        call close_file(f_unit)
#endif
      write(unit, nml=multi_gases_nml)
      call multi_gases_init(ncnst,nwat)

      return
      end subroutine read_namelist_multi_gases_nml


! ----------------------------------------------------------------

! --------------------------------------------------------
      pure real function virq(q)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: variable gas 1+zvir/(1-qc)
!--------------------------------------------
      real, intent(in)           :: q(num_gas)
! Local:
      integer :: n

      virq = vir(sphum)*q(sphum)
      do n=ind_gas,num_gas
         virq = virq+vir(n)*q(sphum+n-1)
      end do
      virq = vir(0)+virq/(1.0-sum(q(sphump1:sphum+num_wat-1)))

      return
      end function
!--------------------------------------------

! --------------------------------------------------------
      pure real function virq_nodq(q)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: variable gas 1+zvir without dividing by 1-qv or 1-qv-qc
!--------------------------------------------
      real, intent(in) :: q(num_gas)
! Local:
      integer :: n

      virq_nodq = vir(0)+vir(sphum)*q(sphum)
      do n=ind_gas,num_gas
         virq_nodq = virq_nodq+vir(n)*q(sphum+n-1)
      end do

      return
      end function
!--------------------------------------------

! --------------------------------------------------------
      pure real function virq_max(q, qmin)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: variable gas 1+zvir using max(qmin,q(sphum))
!--------------------------------------------
      real, intent(in) :: q(num_gas)
      real, intent(in) :: qmin
! Local:
      integer :: n

      virq_max = vir(sphum)*max(qmin,q(sphum))
      do n=ind_gas,num_gas
         virq_max = virq_max+vir(n)*q(sphum+n-1)
      end do
      virq_max = vir(0)+virq_max/(1.0-sum(q(sphump1:sphum+num_wat-1)))


      return
      end function
!--------------------------------------------

! --------------------------------------------------------
      pure real function virq_qpz(q, qpz)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: variable gas 1+zvir/(1.-qpz): qpz in place of qv+qc from q
!--------------------------------------------
      real, intent(in) :: q(num_gas)
      real, intent(in) :: qpz
! Local:
      integer :: n

      virq_qpz = vir(sphum)*q(sphum)
      do n=ind_gas,num_gas
         virq_qpz = virq_qpz+vir(n)*q(sphum+n-1)
      end do
      virq_qpz = vir(0)+virq_qpz/(1.0-qpz)


      return
      end function
!--------------------------------------------

! --------------------------------------------------------
      pure real function virqd(q)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: variable gas 1+zvir/(1-(qv+qc)) (dry)
!--------------------------------------------
      real, intent(in) :: q(num_gas)
! Local:
      integer :: n

      virqd = 0.0
      do n=ind_gas,num_gas
         virqd = virqd+vir(n)*q(sphum+n-1)
      end do
      virqd = vir(0)+virqd/(1.0-sum(q(sphum:sphum+num_wat-1)))

      return
      end function
!--------------------------------------------

! --------------------------------------------------------
      pure real function vicpqd(q)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: variable gas cp (dry)
!--------------------------------------------
      real, intent(in) :: q(num_gas)
! Local:
      integer :: n

      vicpqd = 0.0
      do n=ind_gas,num_gas
         vicpqd = vicpqd+vicp(n)*q(sphum+n-1)
      end do
      vicpqd = vicp(0)+vicpqd/(1.0-sum(q(sphum:sphum+num_wat-1)))

      return
      end function
!--------------------------------------------

! --------------------------------------------------------
      pure real function vicpqd_qpz(q, qpz)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: variable gas cp (dry) with qpz in place of qv+qc from q
!--------------------------------------------
      real, intent(in) :: q(num_gas)
      real, intent(in) :: qpz
! Local:
      integer :: n

      vicpqd_qpz = 0.0
      do n=ind_gas,num_gas
         vicpqd_qpz = vicpqd_qpz+vicp(n)*q(sphum+n-1)
      end do
      vicpqd_qpz = vicp(0)+vicpqd_qpz/(1.0-qpz)

      return
      end function
!--------------------------------------------

! --------------------------------------------------------
      pure real function vicvqd(q)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: variable gas cv (dry)
!--------------------------------------------
      real, intent(in) :: q(num_gas)
! Local:
      integer :: n

      vicvqd = 0.0
      do n=ind_gas,num_gas
         vicvqd = vicvqd+vicv(n)*q(sphum+n-1)
      end do
      vicvqd = vicv(0)+vicvqd/(1.0-sum(q(sphum:sphum+num_wat-1)))

      return
      end function
!--------------------------------------------

! --------------------------------------------------------
      pure real function vicvqd_qpz(q,qpz)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: variable gas cv (dry) with qpz in place of qv+qc from q
!--------------------------------------------
      real, intent(in) :: q(num_gas)
      real, intent(in) :: qpz
! Local:
      integer :: n

      vicvqd_qpz = 0.0
      do n=ind_gas,num_gas
         vicvqd_qpz = vicvqd_qpz+vicv(n)*q(sphum+n-1)
      end do
      vicvqd_qpz = vicv(0)+vicvqd_qpz/(1.0-qpz)

      return
      end function
!--------------------------------------------

end module multi_gases_mod
