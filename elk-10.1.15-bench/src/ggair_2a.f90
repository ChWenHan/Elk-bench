
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_2a
! !INTERFACE:
subroutine ggair_2a(rho,g2rho,gvrho,grho2)
! !USES:
use modmain
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggair\_sp\_2a}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rho(ngtot)
real(8), intent(out) :: g2rho(ngtot),gvrho(ngtot,3),grho2(ngtot)
! local variables
integer ig,ifg,i
! allocatable arrays
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(zfft1(nfgrz),zfft2(nfgrz))
! Fourier transform density to G-space
call rzfftifc(3,ngridg,-1,rho,zfft1)
! grad^2 rho
do ifg=1,nfgrz
  ig=igrzf(ifg)
  if (ig <= ngvc) then
    zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
  else
    zfft2(ifg)=0.d0
  end if
end do
call rzfftifc(3,ngridg,1,g2rho,zfft2)
! grad rho
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    if (ig <= ngvc) then
      zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
    else
      zfft2(ifg)=0.d0
    end if
  end do
  call rzfftifc(3,ngridg,1,gvrho(:,i),zfft2)
end do
! (grad rho)^2
grho2(:)=gvrho(:,1)**2+gvrho(:,2)**2+gvrho(:,3)**2
deallocate(zfft1,zfft2)
end subroutine
!EOC

