
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ggair_3(vx,vc,wx,wc,dtdg2r)
use modmain
implicit none
! arguments
real(8), intent(inout) :: vx(ngtot),vc(ngtot)
real(8), intent(in) :: wx(ngtot),wc(ngtot)
real(8), intent(in) :: dtdg2r(ngtot)
! local variables
integer ig,ifg
! allocatable arrays
real(8), allocatable :: rfir(:)
complex(8), allocatable :: zfft(:)
allocate(rfir(ngtot),zfft(nfgrz))
!------------------!
!     exchange     !
!------------------!
! Fourier transform (wx * dtdg2r) to G-space
rfir(:)=wx(:)*dtdg2r(:)
call rzfftifc(3,ngridg,-1,rfir,zfft)
! grad^2 (wx * dtdg2r)
do ifg=1,nfgrz
  ig=igrzf(ifg)
  zfft(ifg)=-(gc(ig)**2)*zfft(ifg)
end do
call rzfftifc(3,ngridg,1,rfir,zfft)
vx(:)=vx(:)+rfir(:)
!---------------------!
!     correlation     !
!---------------------!
! Fourier transform (wc * dtdg2r) to G-space
rfir(:)=wc(:)*dtdg2r(:)
call rzfftifc(3,ngridg,-1,rfir,zfft)
! grad^2 (wc * dtdg2r)
do ifg=1,nfgrz
  ig=igrzf(ifg)
  zfft(ifg)=-(gc(ig)**2)*zfft(ifg)
end do
call rzfftifc(3,ngridg,1,rfir,zfft)
vc(:)=vc(:)+rfir(:)
deallocate(rfir,zfft)
end subroutine

