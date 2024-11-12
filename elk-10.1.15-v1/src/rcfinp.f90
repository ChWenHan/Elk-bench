
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function rcfinp(rfmt,rfir,cfmt,cfir)
use modmain
use modomp
implicit none
! arguments
real(8), intent(in) :: rfmt(npcmtmax,natmtot),rfir(ngtot)
complex(4), intent(in) :: cfmt(npcmtmax,natmtot),cfir(ngtot)
! local variables
integer is,ias,nthd
! external functions
complex(8), external :: rcfmtinp
! interstitial contribution
rcfinp=sum((cfunir(1:ngtot)*rfir(1:ngtot))*cfir(1:ngtot))
rcfinp=rcfinp*(omega/dble(ngtot))
! muffin-tin contribution
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) REDUCTION(+:rcfinp) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  rcfinp=rcfinp+rcfmtinp(nrcmt(is),nrcmti(is),wrcmt(:,is),rfmt(:,ias), &
   cfmt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end function

