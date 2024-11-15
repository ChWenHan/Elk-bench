
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gencvclmt(nr,nri,ld1,rl,wpr,ld2,crhomt,cvclmt)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: nr(nspecies),nri(nspecies)
integer, intent(in) :: ld1
real(8), intent(in) :: rl(ld1,-lmaxo-1:lmaxo+2,nspecies)
real(8), intent(in) :: wpr(4,ld1,nspecies)
integer, intent(in) :: ld2
complex(4), intent(in) :: crhomt(ld2,natmtot)
complex(4), intent(out) :: cvclmt(ld2,natmtot)
! local variables
integer is,ias,nthd
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call cpotclmt(nr(is),nri(is),ld1,rl(:,:,is),wpr(:,:,is),crhomt(:,ias), &
   cvclmt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

