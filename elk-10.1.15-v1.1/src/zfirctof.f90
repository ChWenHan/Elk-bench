
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfirctof(zfirc,zfir)
use modmain
implicit none
! arguments
complex(8), intent(in) :: zfirc(ngtc)
complex(8), intent(out) :: zfir(ngtot)
! local variables
integer ig
! automatic arrays
complex(8) zfftc(ngtc)
! Fourier transform function on coarse grid to G-space
zfftc(1:ngtc)=zfirc(1:ngtc)
call zfftifc(3,ngdgc,-1,zfftc)
! Fourier transform to fine real-space grid
do ig=1,ngvc
  zfir(igfft(ig))=zfftc(igfc(ig))
end do
zfir(igfft(ngvc+1:ngtot))=0.d0
call zfftifc(3,ngridg,1,zfir)
end subroutine

