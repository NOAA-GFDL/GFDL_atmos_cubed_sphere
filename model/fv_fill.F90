module fv_fill_mod
     implicit none
     public fillz, pfix

contains

 subroutine fillz(im, km, nq, q, dp)
   integer,  intent(in):: im                ! No. of longitudes
   integer,  intent(in):: km                ! No. of levels
   integer,  intent(in):: nq                ! Total number of tracers
   real , intent(in)::  dp(im,km)       ! pressure thickness
   real , intent(inout) :: q(im,km,nq)   ! tracer mixing ratio
! !LOCAL VARIABLES:
   integer i, k, ic
   real  qup, qly, dup

   do ic=1,nq
! Top layer
      do i=1,im
         if( q(i,1,ic) < 0. ) then
             q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
             q(i,1,ic) = 0.
          endif
      enddo

! Interior
      do k=2,km-1
         do i=1,im
         if( q(i,k,ic) < 0. ) then
! Borrow from above
             qup =  max(0., q(i,k-1,ic)*dp(i,k-1) )
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min( 0.5*qly, 0.99*qup )
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
! Borrow from below: q(i,k,ic) is still negative at this stage
             q(i,k+1,ic) = q(i,k+1,ic) - (qly-dup)/dp(i,k+1) 
             q(i,k  ,ic) = 0.
          endif
         enddo
      enddo
 
! Bottom layer
      k = km
      do i=1,im
         if( q(i,k,ic)<0. .and. q(i,k-1,ic)>0.) then
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min(qly, qup)
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
             q(i,k,  ic) = q(i,k,  ic) + dup/dp(i,k  )
#ifdef NON_CONSV_Q
          else
             q(i,km,ic) = 0.
#endif
          endif
      enddo
   enddo
 end subroutine fillz

 subroutine pfix(q, qp, im, ipx, acap, cosp2)

 integer im                  ! Longitudes
 real  acap               ! ???
 real  cosp2              ! ???

 real  q(im)              ! Latitude-level field to adjust
 real  qp(im)             ! Second latitude-level field to adjust (usually pole)

! !OUTPUT PARAMETERS:
 integer ipx                 ! Flag:  0 if Q not change, 1 if changed


! !LOCAL VARIABLES:
 integer i
 real  summ, sump, pmean
 
   summ = 0.
   sump = 0.
   do i=1,im
     summ = summ + q(i)
     sump = sump + qp(i)
   enddo
 
   sump = sump/im
   pmean = (sump*acap + summ*cosp2) / (acap + cosp2*im)
 
   do i=1,im
      q(i) = pmean
      qp(i) = pmean
   enddo
 
   if( qp(1) < 0. ) then
      ipx = 1
   endif

 end subroutine pfix

end module fv_fill_mod
