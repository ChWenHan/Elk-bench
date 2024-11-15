
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_2b
! !INTERFACE:
subroutine ggair_2b(g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
! !USES:
use modmain
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggair\_sp\_2b}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: g2rho(ngtot),gvrho(ngtot,3)
real(8), intent(inout) :: vx(ngtot),vc(ngtot)
real(8), intent(in) :: dxdgr2(ngtot),dcdgr2(ngtot)
! local variables
integer ig,ifg,i
! allocatable arrays
real(8), allocatable :: rfir1(:),rfir2(:)
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(rfir1(ngtot),rfir2(ngtot))
allocate(zfft1(nfgrz),zfft2(nfgrz))
!------------------!
!     exchange     !
!------------------!
! compute grad dxdgr2
call rzfftifc(3,ngridg,-1,dxdgr2,zfft1)
! (grad dxdgr2).(grad rho)
rfir1(:)=0.d0
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call rzfftifc(3,ngridg,1,rfir2,zfft2)
  rfir1(:)=rfir1(:)+rfir2(:)*gvrho(:,i)
end do
vx(:)=vx(:)-2.d0*(rfir1(:)+dxdgr2(:)*g2rho(:))
!---------------------!
!     correlation     !
!---------------------!
! compute grad dcdgr2
call rzfftifc(3,ngridg,-1,dcdgr2,zfft1)
! (grad dcdgr2).(grad rho)
rfir1(:)=0.d0
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call rzfftifc(3,ngridg,1,rfir2,zfft2)
  rfir1(:)=rfir1(:)+rfir2(:)*gvrho(:,i)
end do
vc(:)=vc(:)-2.d0*(rfir1(:)+dcdgr2(:)*g2rho(:))
deallocate(rfir1,rfir2,zfft1,zfft2)
end subroutine
!EOC

