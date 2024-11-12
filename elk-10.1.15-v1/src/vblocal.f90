
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vblocal(vmt,vir,bmt,bir)
use modmain
use modomp
implicit none
! arguments
real(8), intent(out) :: vmt(npcmtmax,natmtot),vir(ngtc)
real(8), intent(out) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtc,ndmag)
! local variables
integer idm,is,ias,nthd
integer nrc,nrci,npc
call holdthd(natmtot+1,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ias,is) &
!$OMP PRIVATE(nrc,nrci,npc,idm) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
! convert muffin-tin Kohn-Sham potential to coarse radial mesh
  call rfmtftoc(nrc,nrci,vsmt(:,ias),vmt(:,ias))
! convert to spherical coordinates
  call rbshtip(nrc,nrci,vmt(:,ias))
! multiply by radial integration weights
  call rfcmtwr(nrc,nrci,wrcmt(:,is),vmt(:,ias))
end do
!$OMP END DO NOWAIT
! multiply interstitial Kohn-Sham potential by characteristic function and
! convert to coarse grid
!$OMP SINGLE
call rfirftoc(vsir,vir)
!$OMP END SINGLE NOWAIT
! repeat for the Kohn-Sham magnetic field
if (spinpol) then
  do idm=1,ndmag
!$OMP DO
    do ias=1,natmtot
      is=idxis(ias)
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
      bmt(1:npc,ias,idm)=bsmt(1:npc,ias,idm)
      call rfcmtwr(nrc,nrci,wrcmt(:,is),bmt(:,ias,idm))
    end do
!$OMP END DO NOWAIT
  end do
! convert Kohn-Sham magnetic field to coarse grid
!$OMP SINGLE
  do idm=1,ndmag
    call rfirftoc(bsir(:,idm),bir(:,idm))
  end do
!$OMP END SINGLE
end if
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

