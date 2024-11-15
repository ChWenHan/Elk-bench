
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine trimrfg(rfir)
use modmain
implicit none
! arguments
real(8), intent(inout) :: rfir(ngtot)
! automatic arrays
complex(8) zfft(nfgrz)
! Fourier transform function to G-space
call rzfftifc(3,ngridg,-1,rfir,zfft)
! zero the components for |G| > 2 gkmax
where(igrzf(1:nfgrz) > ngvc) zfft(1:nfgrz)=0.d0
! Fourier transform back to real-space
call rzfftifc(3,ngridg,1,rfir,zfft)
end subroutine

