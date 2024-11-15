
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_1
! !INTERFACE:
subroutine ggair_1(rho,grho,g2rho,g3rho)
! !USES:
use modmain
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggair\_sp\_1}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rho(ngtot)
real(8), intent(out) :: grho(ngtot),g2rho(ngtot),g3rho(ngtot)
! local variables
integer ig,ifg,i
! allocatable arrays
real(8), allocatable :: gvrho(:,:),rfir(:)
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(gvrho(ngtot,3),rfir(ngtot))
allocate(zfft1(nfgrz),zfft2(nfgrz))
call rzfftifc(3,ngridg,-1,rho,zfft1)
! |grad rho|
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
grho(:)=sqrt(gvrho(:,1)**2+gvrho(:,2)**2+gvrho(:,3)**2)
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
! (grad rho).(grad |grad rho|)
call rzfftifc(3,ngridg,-1,grho,zfft1)
g3rho(:)=0.d0
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    if (ig <= ngvec) then
      zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
    else
      zfft2(ifg)=0.d0
    end if
  end do
  call rzfftifc(3,ngridg,1,rfir,zfft2)
  g3rho(:)=g3rho(:)+gvrho(:,i)*rfir(:)
end do
deallocate(gvrho,rfir,zfft1,zfft2)
end subroutine
!EOC

