
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gencfrc
use modmain
implicit none
! local variables
integer ig,ifg
! automatic arrays
complex(8) zfft(nfgrzc)
do ifg=1,nfgrzc
  ig=igrzfc(ifg)
  if (ig <= ngvc) then
    zfft(ifg)=cfunig(ig)
  else
    zfft(ifg)=0.d0
  end if
end do
! allocate global array
if (allocated(cfrc)) deallocate(cfrc)
allocate(cfrc(ngtc))
! Fourier transform to real-space coarse grid
call rzfftifc(3,ngdgc,1,cfrc,zfft)
end subroutine

