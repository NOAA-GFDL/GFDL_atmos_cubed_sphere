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

!>@brief Molecular diffusion coefficients of
!>       viscosity (mur), conductivity (lam), and diffusivity (d12)
!>       and their efolding time effectiveness
!>@author H.-M. H. Juang, NOAA/NWS/NCEP/EMC

module molecular_diffusion_mod

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

      use constants_mod,      only: rdgas, cp_air
      use fv_mp_mod,          only: is_master
      use mpp_mod,            only: stdlog, input_nml_file
      use fms_mod,            only: check_nml_error, open_namelist_file, close_file
      use fv_grid_utils_mod,  only: g_sum
      use mpp_domains_mod,    only: domain2d
      use fv_arrays_mod,      only: fv_grid_type, fv_grid_bounds_type


#ifdef MULTI_GASES
      use multi_gases_mod,   only: ind_gas, num_gas, virqd, vicpqd, vicvq, vicvqd, virq
#endif


      implicit none

      integer :: ind_gas_str, ind_gas_end
      integer :: md_layers, md_tadj_layers
      real tau_visc, tau_cond, tau_diff
      real md_consv_te
      real md_wait_hr
      real md_wait_sec
      logical md_time
      real, parameter:: amo=15.9994, amo2=31.9999, amo3=47.9982     !g/mol
      real, parameter::              amn2=28.013,  amh2o=18.0154    !g/mol
!< muo3 and muh2o are not precise, correct later
      real, parameter:: muo=3.9e-7, muo2=4.03e-7,  muo3=4.03e-7     !kg/m/s
      real, parameter::             mun2=3.43e-7,  muh2o=3.43e-7    !kg/m/s
!< lao3 is not precise values, but o3_n is very small
      real, parameter:: lao=75.9e-5, lao2=56.e-5,  lao3=36.e-5     !kg/m/s
      real, parameter::              lan2=56.e-5,  lah2o=55.e-5    !kg/m/s
      real, parameter:: cpo=1299.185, cpo2=918.0969, cpo3=820.2391
      real, parameter::               cpn2=1031.108, cph2o=1846.00
      real, parameter:: avgd=6.0221415e23  ! Avogadro constant
      real, parameter:: bz=1.3806505e-23   ! Boltzmann constant J/K
      real, parameter:: a12=9.69e18 ! O-O2 diffusion params
      real, parameter:: s12=0.774
!
      public molecular_diffusion_init, read_namelist_molecular_diffusion_nml
      public molecular_diffusion_coefs, thermosphere_adjustment

      CONTAINS
! --------------------------------------------------------
      subroutine molecular_diffusion_init(ncnst,nwat)
!--------------------------------------------
! molecular diffusion control for each effect from namelist
! Input:
!        ncnst    : number of all prognostic tracers
!        nwat     : number of water and the end location of water
!--------------------------------------------
      integer, intent(in):: ncnst, nwat
!
      md_wait_sec = md_wait_hr * 3600.0
      md_time = .false.

      if( is_master() ) then
        write(*,*) ' molecular_diffusion is on'
        write(*,*) ' molecular_diffusion initial wait seconds ',md_wait_sec
        write(*,*) ' molecular_diffusion number of layers ',md_layers
        write(*,*) ' thermosphere adjustment number of layers ',md_tadj_layers
        write(*,*) ' energy conservation of MD (0:off, 1:on) ',md_consv_te
        write(*,*) ' viscosity    day ',tau_visc,' with effect ',tau_visc
        write(*,*) ' conductivity day ',tau_cond,' with effect ',tau_cond
        write(*,*) ' diffusivity  day ',tau_diff,' with effect ',tau_diff
      endif

#ifdef MULTI_GASES
      ind_gas_str = ind_gas
      ind_gas_end = num_gas
#else
      ind_gas_str = nwat + 1
      ind_gas_end = ncnst
#endif

      return
      end subroutine molecular_diffusion_init

! --------------------------------------------------------
      subroutine read_namelist_molecular_diffusion_nml(nml_filename,ncnst,nwat)

      character(*), intent(IN) :: nml_filename
      integer, intent(IN) :: ncnst, nwat
      integer :: ierr, f_unit, unit, ios

      namelist /molecular_diffusion_nml/ tau_visc, tau_cond, tau_diff, &
                                         md_layers, md_tadj_layers, &
                                         md_wait_hr, md_consv_te

      unit = stdlog()

#ifdef INTERNAL_FILE_NML

      ! Read molecular_diffusion namelist
        read (input_nml_file,molecular_diffusion_nml,iostat=ios)
        ierr = check_nml_error(ios,'molecular_diffusion_nml')

#else
      ! Read molecular_diffusion namelist
        f_unit = open_namelist_file(nml_filename)
        rewind (f_unit)
        read (f_unit,molecular_diffusion_nml,iostat=ios)
        ierr = check_nml_error(ios,'molecular_diffusion_nml')
        call close_file(f_unit)
#endif
      write(unit, nml=molecular_diffusion_nml)
      call molecular_diffusion_init(ncnst,nwat)

      return
      end subroutine read_namelist_molecular_diffusion_nml

! --------------------------------------------------------
      subroutine molecular_diffusion_coefs(dim,plyr,temp,q,mur,lam,d12)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Input: dim : length of field
!        plyr: pressure at given layer (pascal)
!        temp: temperature (K)
!        q   : tracers, needs only gas (kg/kg)
! Ouput: mur : viscosity for momemtum
!        lam : conductivity for thermodynamics
!        d12 : diffusivity for tracers
!--------------------------------------------
      integer, intent(in):: dim
      real, intent(in):: plyr(dim), temp(dim), q(dim,*)
      real, intent(out):: mur(dim), lam(dim), d12(dim)
! Local:
      integer i, n, spo3, spo, spo2
      real   am, fgas, a12bz, avgdbz
      real   qo,  qo2,  qo3,  qn2,  qh2o
      real   o_n, o2_n, o3_n, n2_n, h2o_n
      real   mu, la, t69, rho, cpx, cvx

!constants
      a12bz = a12 * bz
      avgdbz= avgd * bz
      spo3  = ind_gas_str
      spo   = spo3 + 1
      spo2  = spo  + 1
      if( spo.gt.ind_gas_end ) spo=0
      if( spo2.gt.ind_gas_end ) spo2=0

      do n=1,dim
!check
        if(plyr(n).le.0.0) then
           write(*,*) 'ERROR non positive value of plyr ',plyr(n)
           stop
        endif
        if(temp(n).le.0.0) then
           write(*,*) 'ERROR non positive value of temp ',temp(n)
           stop
        endif

        d12(n) = a12bz*temp(n)**s12 * temp(n)/plyr(n)
        if( d12(n).lt.0.0 .and. is_master() ) then
          write(*,*) 'ERROR negative d12 a12bz temp plyr s12 ', &
                      d12(n),a12bz,temp(n),plyr(n),s12
        endif

        fgas = 1.0
        do i=2,ind_gas_str-1
          fgas = fgas - max(0.0,q(n,i))
        enddo
        fgas = max(min(fgas,1.0),1.0E-20)         ! reasonable assured
        fgas = 1./fgas

        qh2o = max(0.0,q(n,1))*fgas
        qo   = 0.0
        qo2  = 0.0
        qo3  = 0.0
        if(spo .ne.0) qo   = max(0.0,q(n,spo ))*fgas
        if(spo2.ne.0) qo2  = max(0.0,q(n,spo2))*fgas
        if(spo3.ne.0) qo3  = max(0.0,q(n,spo3))*fgas
        qn2 = 1.0 - qo - qo2 - qo3 - qh2o

! reasonable values assure
        qo   = max(min(qo ,1.0),0.0)
        qo2  = max(min(qo2,1.0),0.0)
        qo3  = max(min(qo3,1.0),0.0)
        qn2  = max(min(qn2,1.0),0.0)
        qh2o = max(min(qh2o,1.0),0.0)
        cpx  = ( cpo*qo + cpo2*qo2 + cpo3*qo3 + cpn2*qn2 + cph2o*qh2o ) / &
               (     qo +      qo2 +      qo3 +      qn2 +       qh2o )
        cvx  = cpx - rdgas

        am = qo/amo + qo2/amo2 + qo3/amo3 + qn2/amn2 +qh2o/amh2o
        am = 1.0 / am    ! g/mol
        o_n   = qo   * am / amo
        o2_n  = qo2  * am / amo2
        o3_n  = qo3  * am / amo3
        n2_n  = qn2  * am / amn2
        h2o_n = qh2o * am / amh2o
        mu = o_n*muo + o2_n*muo2 + o3_n*muo3 + n2_n*mun2 + h2o_n*muh2o
        la = o_n*lao + o2_n*lao2 + o3_n*lao3 + n2_n*lan2 + h2o_n*lah2o

        t69 = temp(n) ** 0.69
        rho = 1e-3 * am * plyr(n)/temp(n) / avgdbz    ! km/m^3

        mur(n) = mu * t69 / rho
        lam(n) = la * t69 / rho / cvx
!       lam(n) = la * t69 / rho / cpx

!reasonable assured
         mur(n) = min(mur(n),1.0e15)    ! viscosity
         lam(n) = min(lam(n),1.0e15)    ! conductivity
         d12(n) = min(d12(n),1.0e15)    ! diffusivity

      enddo

      return
      end subroutine molecular_diffusion_coefs
      subroutine thermosphere_adjustment(domain,gridstruct,npz,bd,ng,pt)
      type(fv_grid_bounds_type), intent(IN) :: bd
      real ,      intent(INOUT) ::    pt(bd%is:bd%ie  ,bd%js:bd%je, npz)
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain
      integer, intent(IN) :: ng, npz
      real       :: correct, aave_t1, ttsave, tdsave
      integer :: k, j, i
      real :: t(bd%is:bd%ie,bd%js:bd%je)
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




    k=md_tadj_layers
    do j=js,je
       do i=is,ie
          t(i,j) = pt(i,j,k)
       enddo
    enddo
    aave_t1 = g_sum(domain,t,is,ie,js,je,ng,gridstruct%area_64, 1)
    k=md_tadj_layers-1
    do j=js,je
       do i=is,ie
          t(i,j) = pt(i,j,k)
       enddo
    enddo
    ttsave = g_sum(domain,t,is,ie,js,je,ng,gridstruct%area_64, 1)
    tdsave = ttsave - aave_t1
!   if( is_master() ) write(*,*) ' k t1 ts ',k,aave_t1,ttsave
!   correct_sum = 0.0
    do k=md_tadj_layers-2,1,-1
       do j=js,je
          do i=is,ie
             t(i,j) = pt(i,j,k)
          enddo
       enddo
       aave_t1 = g_sum(domain,t,is,ie,js,je,ng,gridstruct%area_64, 1)
!      if( is_master() ) write(*,*) ' k ts t1 ',k,ttsave,aave_t1
       if( aave_t1 .gt. ttsave ) then
          tdsave = min(tdsave,aave_t1-ttsave)
          ttsave = ttsave + tdsave
       else
          tdsave = 0.0
       endif
       correct = ttsave - aave_t1
!      if( is_master() ) write(*,*) ' k t1 ts correct ',k,aave_t1,ttsave,correct
       if( correct .ne. 0.0 ) then
          do j=js,je
             do i=is,ie
                pt(i,j,k) = t(i,j)  + correct
             enddo
          enddo
       endif
!      correct_sum = correct_sum + correct
    enddo
! adjust for energy conserving by evenly distribute total correction
!   if( correct_sum .ne. 0.0 ) then
!      correct = - correct_sum / (md_tadj_layers + 10)
!      if( is_master() ) write(*,*) ' sum=',correct_sum,' correct=',correct
!      do k=1,md_tadj_layers+10
!         do j=js,je
!            do i=is,ie
!               pt(i,j,k) = t(i,j)  + correct
!            enddo
!         enddo
!      enddo
!   endif
  return
  end subroutine thermosphere_adjustment

end module molecular_diffusion_mod
