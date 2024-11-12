
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine genbs
use modmain
use modomp
implicit none
! local variables
integer idm,is,ia,ias
integer nrc,nrci,npc,nthd
real(8) cb,t1
if (.not.spinpol) return
! coupling constant of the external field (g_e/4c)
cb=gfacte/(4.d0*solsc)
!------------------------------------!
!     muffin-tin Kohn-Sham field     !
!------------------------------------!
call holdthd(natmtot+ndmag,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ias,is,ia,nrc) &
!$OMP PRIVATE(nrci,npc,idm,t1) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
! exchange-correlation magnetic field in spherical coordinates
  do idm=1,ndmag
    call rfmtftoc(nrc,nrci,bxcmt(:,ias,idm),bsmt(:,ias,idm))
    call rbshtip(nrc,nrci,bsmt(:,ias,idm))
  end do
! add the external magnetic field
  t1=cb*(bfcmt(3,ia,is)+bfieldc(3))
  bsmt(1:npc,ias,ndmag)=bsmt(1:npc,ias,ndmag)+t1
  if (ncmag) then
    do idm=1,2
      t1=cb*(bfcmt(idm,ia,is)+bfieldc(idm))
      bsmt(1:npc,ias,idm)=bsmt(1:npc,ias,idm)+t1
    end do
  end if
end do
!$OMP END DO NOWAIT
!-----------------------------------------------!
!     interstitial Kohn-Sham magnetic field     !
!-----------------------------------------------!
!$OMP DO
do idm=1,ndmag
  if (ncmag) then
    t1=cb*bfieldc(idm)
  else
    t1=cb*bfieldc(3)
  end if
  bsir(1:ngtot,idm)=bxcir(1:ngtot,idm)+t1
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
! add the magnetic dipole field if required
if (tbdip) call bdipole
! multiply interstitial part by characteristic function and store on coarse grid
do idm=1,ndmag
  call rfirftoc(bsir(:,idm),bsirc(:,idm))
end do
end subroutine

