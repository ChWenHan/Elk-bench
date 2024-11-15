
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rfirctof
! !INTERFACE:
subroutine rfirctof(rfirc,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   rfirc : real input function on coarse grid (in,real(ngtc))
!   rfir  : real output function on fine grid (out,real(ngtot))
! !DESCRIPTION:
!   Converts a real function on a coarse grid given by sizes {\tt ngdgc} to a
!   function on a fine grid given by {\tt ngridg}. This is done by first
!   Fourier transforming {\tt rfirc} to $G$-space, zeroing the missing values
!   and then transforming back to {\tt rfir}.
!
! !REVISION HISTORY:
!   Created March 2020 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rfirc(ngtc)
real(8), intent(out) :: rfir(ngtot)
! local variables
integer ig,ifg
! automatic arrays
complex(8) zfftc(ngtc)
! allocatable arrays
complex(8), allocatable :: zfft(:)
! Fourier transform function on coarse grid to G-space
zfftc(1:ngtc)=rfirc(1:ngtc)
call zfftifc(3,ngdgc,-1,zfftc)
! Fourier transform to fine real-space grid
allocate(zfft(nfgrz))
do ifg=1,nfgrz
  ig=igrzf(ifg)
  if (ig <= ngvc) then
    zfft(ifg)=zfftc(igfc(ig))
  else
    zfft(ifg)=0.d0
  end if
end do
call rzfftifc(3,ngridg,1,rfir,zfft)
deallocate(zfft)
end subroutine
!EOC

