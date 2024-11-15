
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine cfmtconj(nr,nri,np,cfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri,np
complex(4), intent(inout) :: cfmt(np)
! local variables
integer i
! automatic arrays
complex(4) cfmt1(np)
cfmt1(1:np)=cfmt(1:np)
call cflmnconj(lmaxi,nri,lmmaxi,cfmt1,cfmt)
i=lmmaxi*nri+1
call cflmnconj(lmaxo,nr-nri,lmmaxo,cfmt1(i),cfmt(i))
return

contains

!BOP
! !ROUTINE: cflmnconj
! !INTERFACE:
pure subroutine cflmnconj(lmax,n,ld,cflm1,cflm2)
! !INPUT/OUTPUT PARAMETERS:
!   lmax  : maximum angular momentum (in,integer)
!   n     : number of functions to conjugate (in,integer)
!   ld    : leading dimension (in,integer)
!   cflm1 : coefficients of input complex spherical harmonic expansion
!           (in,complex((lmax+1)**2)))
!   cflm2 : coefficients of output complex spherical harmonic expansion
!           (out,complex((lmax+1)**2)))
! !DESCRIPTION:
!   Returns the complex conjugate of a function expanded in spherical harmonics.
!   In other words, given the input function coefficients $c_{lm}$, the routine
!   returns  $c'_{lm}=(-1)^m c^*_{l-m}$ so that
!   $$ \sum_{lm}c'_{lm}Y_{lm}(\theta,\phi)=\left(\sum_{lm}c_{lm}Y_{lm}
!    (\theta,\phi)\right)^* $$
!   for all $(\theta,\phi)$.
!
! !REVISION HISTORY:
!   Created April 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax,n,ld
complex(4), intent(in) :: cflm1(ld,n)
complex(4), intent(out) :: cflm2(ld,n)
! local variables
integer l,m,lm1,lm2
do l=0,lmax
  lm1=l**2
  lm2=(l+1)**2+1
  do m=-l,-1
    lm1=lm1+1
    lm2=lm2-1
    if (mod(m,2) == 0) then
      cflm2(lm1,1:n)=conjg(cflm1(lm2,1:n))
      cflm2(lm2,1:n)=conjg(cflm1(lm1,1:n))
    else
      cflm2(lm1,1:n)=-conjg(cflm1(lm2,1:n))
      cflm2(lm2,1:n)=-conjg(cflm1(lm1,1:n))
    end if
  end do
! m=0 case
  lm1=lm1+1
  cflm2(lm1,1:n)=conjg(cflm1(lm1,1:n))
end do
end subroutine
!EOC

end subroutine

