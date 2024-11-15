
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfirftoc(zfir,zfirc)
use modmain
implicit none
! arguments
complex(8), intent(in) :: zfir(ngtot)
complex(8), intent(out) :: zfirc(ngtc)
! local variables
integer ig
! allocatable arrays
complex(8), allocatable :: zfft(:)
allocate(zfft(ngtot))
! multiply by characteristic function and Fourier transform on fine grid
zfft(1:ngtot)=zfir(1:ngtot)*cfunir(1:ngtot)
call zfftifc(3,ngridg,-1,zfft)
zfirc(1:ngtc)=0.d0
do ig=1,ngvc
  zfirc(igfc(ig))=zfft(igfft(ig))
end do
! Fourier transform to real-space coarse grid
call zfftifc(3,ngdgc,1,zfirc)
deallocate(zfft)
end subroutine

