
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genylmg
! !INTERFACE:
subroutine genylmg
! !USES:
use modmain
! !DESCRIPTION:
!   Generates a set of spherical harmonics, $4\pi(-i)^l Y_{lm}(\hat{\bf G})$,
!   with angular momenta up to {\tt lmaxo} for the set of ${\bf G}$-vectors.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ig
! allocate global G-vector spherical harmonic array
if (allocated(ylmg)) deallocate(ylmg)
allocate(ylmg(lmmaxo,ngvec))
do ig=1,ngvec
  call genylmv(.true.,lmaxo,vgc(:,ig),ylmg(:,ig))
end do
end subroutine
!EOC

