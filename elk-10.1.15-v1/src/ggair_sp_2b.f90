
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_sp_2b
! !INTERFACE:
subroutine ggair_sp_2b(g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2,dxdgd2, &
 dxdgud,dcdgu2,dcdgd2,dcdgud)
! !USES:
use modmain
! !DESCRIPTION:
!   Post processing step of interstitial gradients for GGA type 2. See routine
!   {\tt ggamt\_sp\_2a} for full details.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
real(8), intent(in) :: g2up(ngtot),g2dn(ngtot)
real(8), intent(in) :: gvup(ngtot,3),gvdn(ngtot,3)
real(8), intent(inout) :: vxup(ngtot),vxdn(ngtot)
real(8), intent(inout) :: vcup(ngtot),vcdn(ngtot)
real(8), intent(in) :: dxdgu2(ngtot),dxdgd2(ngtot),dxdgud(ngtot)
real(8), intent(in) :: dcdgu2(ngtot),dcdgd2(ngtot),dcdgud(ngtot)
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
! compute grad dxdgu2
call rzfftifc(3,ngridg,-1,dxdgu2,zfft1)
! (grad dxdgu2).(grad rhoup)
rfir1(:)=0.d0
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call rzfftifc(3,ngridg,1,rfir2,zfft2)
  rfir1(:)=rfir1(:)+rfir2(:)*gvup(:,i)
end do
vxup(:)=vxup(:)-2.d0*(rfir1(:)+dxdgu2(:)*g2up(:))-dxdgud(:)*g2dn(:)
! compute grad dxdgd2
call rzfftifc(3,ngridg,-1,dxdgd2,zfft1)
! (grad dxdgd2).(grad rhodn)
rfir1(:)=0.d0
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call rzfftifc(3,ngridg,1,rfir2,zfft2)
  rfir1(:)=rfir1(:)+rfir2(:)*gvdn(:,i)
end do
vxdn(:)=vxdn(:)-2.d0*(rfir1(:)+dxdgd2(:)*g2dn(:))-dxdgud(:)*g2up(:)
! compute grad dxdgud
call rzfftifc(3,ngridg,-1,dxdgud,zfft1)
! (grad dxdgud).(grad rhodn) and (grad dxdgud).(grad rhoup)
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call rzfftifc(3,ngridg,1,rfir2,zfft2)
  vxup(:)=vxup(:)-rfir2(:)*gvdn(:,i)
  vxdn(:)=vxdn(:)-rfir2(:)*gvup(:,i)
end do
!---------------------!
!     correlation     !
!---------------------!
! compute grad dcdgu2
call rzfftifc(3,ngridg,-1,dcdgu2,zfft1)
! (grad dcdgu2).(grad rhoup)
rfir1(:)=0.d0
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call rzfftifc(3,ngridg,1,rfir2,zfft2)
  rfir1(:)=rfir1(:)+rfir2(:)*gvup(:,i)
end do
vcup(:)=vcup(:)-2.d0*(rfir1(:)+dcdgu2(:)*g2up(:))-dcdgud(:)*g2dn(:)
! compute grad dcdgd2
call rzfftifc(3,ngridg,-1,dcdgd2,zfft1)
! (grad dcdgd2).(grad rhodn)
rfir1(:)=0.d0
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call rzfftifc(3,ngridg,1,rfir2,zfft2)
  rfir1(:)=rfir1(:)+rfir2(:)*gvdn(:,i)
end do
vcdn(:)=vcdn(:)-2.d0*(rfir1(:)+dcdgd2(:)*g2dn(:))-dcdgud(:)*g2up(:)
! compute grad dcdgud
call rzfftifc(3,ngridg,-1,dcdgud,zfft1)
! (grad dcdgud).(grad rhodn) and (grad dcdgud).(grad rhoup)
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call rzfftifc(3,ngridg,1,rfir2,zfft2)
  vcup(:)=vcup(:)-rfir2(:)*gvdn(:,i)
  vcdn(:)=vcdn(:)-rfir2(:)*gvup(:,i)
end do
deallocate(rfir1,rfir2,zfft1,zfft2)
end subroutine
!EOC

