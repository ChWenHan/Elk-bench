
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_sp_2a
! !INTERFACE:
subroutine ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the interstitial gradients $\nabla^2\rho^{\uparrow}$,
!   $\nabla^2\rho^{\downarrow}$, $\nabla\rho^{\uparrow}$,
!   $\nabla\rho^{\downarrow}$, $(\nabla\rho^{\uparrow})^2$,
!   $(\nabla\rho^{\downarrow})^2$ and
!   $\nabla\rho^{\uparrow}\cdot\nabla\rho^{\downarrow}$. These are used for GGA
!   functionals of type 2 and meta-GGA. See {\tt ggamt\_sp\_2a} for details.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rhoup(ngtot),rhodn(ngtot)
real(8), intent(out) :: g2up(ngtot),g2dn(ngtot)
real(8), intent(out) :: gvup(ngtot,3),gvdn(ngtot,3)
real(8), intent(out) :: gup2(ngtot),gdn2(ngtot),gupdn(ngtot)
! local variables
integer ig,ifg,i
! allocatable arrays
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(zfft1(nfgrz),zfft2(nfgrz))
!----------------!
!     rho up     !
!----------------!
call rzfftifc(3,ngridg,-1,rhoup,zfft1)
! compute grad^2 rhoup
do ifg=1,nfgrz
  ig=igrzf(ifg)
  if (ig <= ngvc) then
    zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
  else
    zfft2(ifg)=0.d0
  end if
end do
call rzfftifc(3,ngridg,1,g2up,zfft2)
! grad rhoup
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    if (ig <= ngvc) then
      zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
    else
      zfft2(ifg)=0.d0
    end if
  end do
  call rzfftifc(3,ngridg,1,gvup(:,i),zfft2)
end do
! (grad rhoup)^2
gup2(:)=gvup(:,1)**2+gvup(:,2)**2+gvup(:,3)**2
!------------------!
!     rho down     !
!------------------!
call rzfftifc(3,ngridg,-1,rhodn,zfft1)
! compute grad^2 rhodn
do ifg=1,nfgrz
  ig=igrzf(ifg)
  if (ig <= ngvc) then
    zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
  else
    zfft2(ifg)=0.d0
  end if
end do
call rzfftifc(3,ngridg,1,g2dn,zfft2)
! grad rhodn
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    if (ig <= ngvc) then
      zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
    else
      zfft2(ifg)=0.d0
    end if
  end do
  call rzfftifc(3,ngridg,1,gvdn(:,i),zfft2)
end do
! (grad rhodn)^2
gdn2(:)=gvdn(:,1)**2+gvdn(:,2)**2+gvdn(:,3)**2
! (grad rhoup).(grad rhodn)
gupdn(:)=gvup(:,1)*gvdn(:,1)+gvup(:,2)*gvdn(:,2)+gvup(:,3)*gvdn(:,3)
deallocate(zfft1,zfft2)
end subroutine
!EOC

