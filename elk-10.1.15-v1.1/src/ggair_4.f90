
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ggair_4(gvrho,vx,vc,wx,wc,dtdr,dtdgr2)
use modmain
implicit none
! arguments
real(8), intent(in) :: gvrho(ngtot,3)
real(8), intent(inout) :: vx(ngtot),vc(ngtot)
real(8), intent(in) :: wx(ngtot),wc(ngtot)
real(8), intent(in) :: dtdr(ngtot),dtdgr2(ngtot)
! local variables
integer ig,ifg,i
! allocatable arrays
real(8), allocatable :: rfir1(:),rfir2(:)
complex(8), allocatable :: zfft(:)
allocate(rfir1(ngtot),rfir2(ngtot),zfft(nfgrz))
!------------------!
!     exchange     !
!------------------!
vx(:)=vx(:)+wx(:)*dtdr(:)
rfir1(:)=wx(:)*dtdgr2(:)
do i=1,3
  rfir2(:)=rfir1(:)*gvrho(:,i)
  call rzfftifc(3,ngridg,-1,rfir2,zfft)
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    zfft(ifg)=vgc(i,ig)*cmplx(-aimag(zfft(ifg)),dble(zfft(ifg)),8)
  end do
  call rzfftifc(3,ngridg,1,rfir2,zfft)
  vx(:)=vx(:)-2.d0*rfir2(:)
end do
!---------------------!
!     correlation     !
!---------------------!
vc(:)=vc(:)+wc(:)*dtdr(:)
rfir1(:)=wc(:)*dtdgr2(:)
do i=1,3
  rfir2(:)=rfir1(:)*gvrho(:,i)
  call rzfftifc(3,ngridg,-1,rfir2,zfft)
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    zfft(ifg)=vgc(i,ig)*cmplx(-aimag(zfft(ifg)),dble(zfft(ifg)),8)
  end do
  call rzfftifc(3,ngridg,1,rfir2,zfft)
  vc(:)=vc(:)-2.d0*rfir2(:)
end do
deallocate(rfir1,rfir2,zfft)
end subroutine

