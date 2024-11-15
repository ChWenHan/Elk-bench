
! Copyright (C) 2017 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gwlocal(vmt,vir,bmt,bir)
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
  npc=npcmt(is)
! convert exchange-correlation potential to a coarse radial mesh
  call rfmtftoc(nrc,nrci,vxcmt(:,ias),vmt(:,ias))
! negate because V_xc should be removed from the self-energy
  vmt(1:npc,ias)=-vmt(1:npc,ias)
! convert to spherical coordinates
  call rbshtip(nrc,nrci,vmt(:,ias))
! multiply by radial integration weights
  call rfcmtwr(nrc,nrci,wrcmt(:,is),vmt(:,ias))
end do
!$OMP END DO NOWAIT
! multiply -V_xc by the characteristic function and convert to coarse grid
!$OMP SINGLE
call rfirftoc(vxcir,vir)
vir(1:ngtc)=-vir(1:ngtc)
!$OMP END SINGLE NOWAIT
! do the same for B_xc in the spin-polarised case
if (spinpol) then
  do idm=1,ndmag
!$OMP DO
    do ias=1,natmtot
      is=idxis(ias)
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
      call rfmtftoc(nrc,nrci,bxcmt(:,ias,idm),bmt(:,ias,idm))
      bmt(1:npc,ias,idm)=-bmt(1:npc,ias,idm)
      call rbshtip(nrc,nrci,bmt(:,ias,idm))
      call rfcmtwr(nrc,nrci,wrcmt(:,is),bmt(:,ias,idm))
    end do
!$OMP END DO NOWAIT
  end do
!$OMP SINGLE
  do idm=1,ndmag
    call rfirftoc(bxcir(:,idm),bir(:,idm))
    bir(1:ngtc,idm)=-bir(1:ngtc,idm)
  end do
!$OMP END SINGLE
end if
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

