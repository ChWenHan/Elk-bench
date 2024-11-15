
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gencfrm(wfmt11,wfmt12,wfir11,wfir12,wfmt21,wfmt22,wfir21,wfir22, &
 crhomt,crhoir,cmagmt,cmagir)
use modmain
use modomp
implicit none
! arguments
complex(4), intent(in) ::  wfmt11(npcmtmax,natmtot),wfmt12(npcmtmax,natmtot)
complex(4), intent(in) ::  wfir11(ngtot),wfir12(ngtot)
complex(4), intent(in) ::  wfmt21(npcmtmax,natmtot),wfmt22(npcmtmax,natmtot)
complex(4), intent(in) ::  wfir21(ngtot),wfir22(ngtot)
complex(4), intent(out) :: crhomt(npcmtmax,natmtot),crhoir(ngtot)
complex(4), intent(out) :: cmagmt(npcmtmax,natmtot,ndmag),cmagir(ngtot,ndmag)
! local variables
integer ld,is,ias,nthd
ld=npcmtmax*natmtot
call holdthd(natmtot+1,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ias,is) &
!$OMP NUM_THREADS(nthd)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call gencrm(npcmt(is),wfmt11(:,ias),wfmt12(:,ias),wfmt21(:,ias), &
   wfmt22(:,ias),crhomt(:,ias),ld,cmagmt(:,ias,1))
end do
!$OMP END DO NOWAIT
!---------------------------!
!     interstitial part     !
!---------------------------!
!$OMP SINGLE
call gencrm(ngtot,wfir11,wfir12,wfir21,wfir22,crhoir,ngtot,cmagir)
!$OMP END SINGLE
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

