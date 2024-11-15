
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfirftoc(rfir,rfirc)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfir(ngtot)
real(8), intent(out) :: rfirc(ngtc)
! local variables
integer ig,ifg
! allocatable arrays
complex(8), allocatable :: zfft(:)
! automatic arrays
complex(8) zfftc(nfgrzc)
allocate(zfft(ngtot))
! multiply by characteristic function and Fourier transform on fine grid
zfft(1:ngtot)=rfir(1:ngtot)*cfunir(1:ngtot)
call zfftifc(3,ngridg,-1,zfft)
! Fourier transform to coarse real-space grid
do ifg=1,nfgrzc
  ig=igrzfc(ifg)
  if (ig <= ngvc) then
    zfftc(ifg)=zfft(igfft(ig))
  else
    zfftc(ifg)=0.d0
  end if
end do
! Fourier transform to real-space coarse grid
call rzfftifc(3,ngdgc,1,rfirc,zfftc)
deallocate(zfft)
end subroutine

