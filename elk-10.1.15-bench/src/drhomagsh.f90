
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine drhomagsh
use modmain
use modphonon
use modomp
implicit none
! local variables
integer idm,is,ias,nthd
call holdthd(2*natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is,idm) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
! convert the density derivative to spherical harmonics
  call zfshtip(nrcmt(is),nrcmti(is),drhomt(:,ias))
end do
!$OMP END DO NOWAIT
do idm=1,ndmag
!$OMP DO
  do ias=1,natmtot
    is=idxis(ias)
! convert the magnetisation derivative to spherical harmonics
    call zfshtip(nrcmt(is),nrcmti(is),dmagmt(:,ias,idm))
  end do
!$OMP END DO
end do
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

