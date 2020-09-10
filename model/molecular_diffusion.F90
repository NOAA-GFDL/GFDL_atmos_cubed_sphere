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

      use constants_mod,     only: rdgas, cp_air
      use     fv_mp_mod,     only: is_master
#ifdef MULTI_GASES
      use multi_gases_mod,   only: ind_gas, num_gas
#endif


      implicit none

      integer :: ind_gas_str, ind_gas_end
      integer :: md_layers
      logical md_init_wait
      real atau_visc, atau_cond, atau_diff    
      real md_wait_sec
      real, parameter:: amo=15.9994, amo2=31.9999, amo3=47.9982     !g/mol
      real, parameter::              amn2=28.013,  amh2o=18.0154    !g/mol
!hmhj muo3 and muh2o are not precise, correct later
      real, parameter:: muo=3.9e-7, muo2=4.03e-7,  muo3=4.03e-7     !kg/m/s
      real, parameter::             mun2=3.43e-7,  muh2o=3.43e-7    !kg/m/s
! hmhj lao3 is not precise values, but o3_n is very small
      real, parameter:: lao=75.9e-5, lao2=56.e-5,  lao3=36.e-5     !kg/m/s
      real, parameter::              lan2=56.e-5,  lah2o=55.e-5    !kg/m/s
      real, parameter:: cpo=1299.185, cpo2=918.0969, cpo3=820.2391
      real, parameter::               cpn2=1031.108, cph2o=1846.00
      real, parameter:: avgd=6.0221415e23  ! Avogadro constant
      real, parameter:: bz=1.3806505e-23   ! Boltzmann constant J/K
      real, parameter:: a12=9.69e18 ! O-O2 diffusion params
      real, parameter:: s12=0.774
!
      public molecular_diffusion_init
      public molecular_diffusion_coefs

      CONTAINS
! --------------------------------------------------------
      subroutine molecular_diffusion_init(tau_visc,tau_cond,tau_diff, &
                 md_n_layer,md_wait_hr,ncnst,nwat)
!--------------------------------------------
! molecular diffusion control for each effect
! Input: tau_visc : viscosity effect weighting
!        tau_cond : conductivity effect weighting
!        tau_diff : diffusivity effect weighting
!        md_n_layer  : how many layers are applied md, counted from top
!        md_wait_hr  : 0 for no wait to start otherwise by hour
!        ncnst    : number of all prognostic tracers
!        nwat     : number of water and the end location of water
!--------------------------------------------
      real, intent(in):: tau_visc, tau_cond, tau_diff, md_wait_hr
      integer, intent(in):: md_n_layer, ncnst, nwat
!
      atau_visc = abs(tau_visc)
      atau_cond = abs(tau_cond)
      atau_diff = abs(tau_diff)
      md_layers   = md_n_layer
!
      if( md_wait_hr.gt.0.0 ) then
        md_init_wait= .true.
        md_wait_sec = md_wait_hr * 3600.0
      else
        md_init_wait= .false.
        md_wait_sec = 0.0
      endif
      
      if( is_master() ) then
        write(*,*) ' molecular_diffusion is on'
        write(*,*) ' molecular_diffusion initial wait is ',md_init_wait
        write(*,*) ' molecular_diffusion initial wait seconds ',md_wait_sec
        write(*,*) ' molecular_diffusion number of layers ',md_layers
        write(*,*) ' viscosity    day ',tau_visc,' with effect ',atau_visc
        write(*,*) ' conductivity day ',tau_cond,' with effect ',atau_cond
        write(*,*) ' diffusivity  day ',tau_diff,' with effect ',atau_diff
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
      integer i, n, spfo3, spfo, spfo2
      real   am, fgas, a12bz, avgdbz
      real   qo,  qo2,  qo3,  qn2,  qh2o
      real   o_n, o2_n, o3_n, n2_n, h2o_n
      real   mu, la, t69, rho, cpx

!constants
      a12bz = a12 * bz
      avgdbz= avgd * bz
      spfo3 = ind_gas_str
      spfo  = spfo3 + 1
      spfo2 = spfo  + 1
      if( spfo.gt.ind_gas_end ) spfo=0
      if( spfo2.gt.ind_gas_end ) spfo2=0

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
        fgas = max(min(fgas,1.0),1.0E-20)	! reasonable assured
        fgas = 1./fgas

        qh2o = max(0.0,q(n,1))*fgas
        qo   = 0.0
        qo2  = 0.0
        qo3  = 0.0
        if(spfo .ne.0) qo   = max(0.0,q(n,spfo ))*fgas
        if(spfo2.ne.0) qo2  = max(0.0,q(n,spfo2))*fgas
        if(spfo3.ne.0) qo3  = max(0.0,q(n,spfo3))*fgas
        qn2 = 1.0 - qo - qo2 - qo3 - qh2o

! reasonable values assure
        qo   = max(min(qo ,1.0),0.0)
        qo2  = max(min(qo2,1.0),0.0)
        qo3  = max(min(qo3,1.0),0.0)
        qn2  = max(min(qn2,1.0),0.0)
        qh2o = max(min(qh2o,1.0),0.0)
        cpx  = ( cpo*qo + cpo2*qo2 + cpn2*qn2 + cph2o*qh2o ) / &
               (     qo +      qo2 +      qn2 +       qh2o )

        am = qo/amo + qo2/amo2 + qo3/amo3 + qn2/amn2 +qh2o/amh2o
        am = 1.0 / am				! g/mol
        o_n   = qo   * am / amo
        o2_n  = qo2  * am / amo2
        o3_n  = qo3  * am / amo3
        n2_n  = qn2  * am / amn2
        h2o_n = qh2o * am / amh2o
        mu = o_n*muo + o2_n*muo2 + o3_n*muo3 + n2_n*mun2 + h2o_n*muh2o
        la = o_n*lao + o2_n*lao2 + o3_n*lao3 + n2_n*lan2 + h2o_n*lah2o

        t69 = temp(n) ** 0.69
        rho = 1e-3 * am * plyr(n)/temp(n) / avgdbz	! km/m^3
        
        mur(n) = mu * t69 / rho
        lam(n) = la * t69 / rho / cpx

!reasonable assured
         mur(n) = min(mur(n),1.0e15)	! viscosity
         lam(n) = min(lam(n),1.0e15)	! conductivity
         d12(n) = min(d12(n),1.0e15)	! diffusivity

      enddo

      return
      end subroutine molecular_diffusion_coefs
!--------------------------------------------

end module molecular_diffusion_mod
