
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genvsig
! !INTERFACE:
subroutine genvsig
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the Fourier transform of the Kohn-Sham effective potential in the
!   interstitial region.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! automatic arrays
complex(8) zfft(ngtc)
! Fourier transform intersitial potential to G-space
zfft(1:ngtc)=vsirc(1:ngtc)
call zfftifc(3,ngdgc,-1,zfft)
! store in global array
vsig(1:ngvc)=zfft(igfc(1:ngvc))
end subroutine
!EOC

