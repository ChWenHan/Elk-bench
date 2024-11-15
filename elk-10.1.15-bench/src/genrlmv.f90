
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: genrlmv
! !INTERFACE:
subroutine genrlmv(lmax,v,rlm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   v    : input vector (in,real(3))
!   rlm  : array of real spherical harmonics (out,real((lmax+1)**2))
! !DESCRIPTION:
!   Generates a sequence of real spherical harmonics evaluated at angles
!   $(\theta,\phi)$ for $0<l<l_{\rm max}$. The values are returned in a packed
!   array {\tt rlm} indexed with $j=l(l+1)+m+1$. Real spherical harmonics are
!   defined by
!   $$ R_{lm}(\theta,\phi)=\begin{cases}
!     \sqrt{2}\,{\rm Re}\{Y_{lm}(\theta,\phi)\} & m>0 \\
!     \sqrt{2}\,{\rm Im}\{Y_{lm}(\theta,\phi)\} & m<0 \\
!     {\rm Re}\{Y_{lm}(\theta,\phi)\} & m=0
!    \end{cases} $$
!   where $Y_{lm}$ are the complex spherical harmonics. These functions are
!   orthonormal and complete and may be used for expanding real-valued functions
!   on the sphere. This routine is numerically stable and accurate to near
!   machine precision for $l\le 50$. See routine {\tt genylmv}.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
real(8), intent(in) :: v(3)
real(8), intent(out) :: rlm(*)
! local variables
integer l,lm,n
real(8), parameter :: sqtwo=1.4142135623730950488d0
! automatic arrays
complex(8) ylm((lmax+1)**2)
if ((lmax < 0).or.(lmax > 50)) then
  write(*,*)
  write(*,'("Error(genrlmv): lmax out of range : ",I8)') lmax
  write(*,*)
  stop
end if
! generate complex spherical harmonics
call genylmv(.false.,lmax,v,ylm)
! convert to real spherical harmonics
rlm(1)=dble(ylm(1))
do l=1,lmax
  n=l-1
  lm=l**2+1
  rlm(lm:lm+n)=sqtwo*aimag(ylm(lm:lm+n))
  lm=lm+l
  rlm(lm)=dble(ylm(lm))
  lm=lm+1
  rlm(lm:lm+n)=sqtwo*dble(ylm(lm:lm+n))
end do
end subroutine
!EOC

