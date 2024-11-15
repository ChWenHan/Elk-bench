
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendvsig
use modmain
use modphonon
implicit none
! allocatable arrays
complex(8), allocatable :: zfft(:)
allocate(zfft(ngtot))
zfft(1:ngtot)=dvsir(1:ngtot)*cfunir(1:ngtot)
zfft(1:ngtot)=zfft(1:ngtot)+vsir(1:ngtot)*dcfunir(1:ngtot)
! Fourier transform to G+q-space
call zfftifc(3,ngridg,-1,zfft)
! store in global array
dvsig(1:ngvc)=zfft(igfft(1:ngvc))
deallocate(zfft)
end subroutine

