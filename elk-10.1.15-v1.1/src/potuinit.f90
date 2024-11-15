
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potuinit
use modmain
use modulr
use modrandom
use modomp
implicit none
! local variables
integer ifq,idm,is,ias
integer nrc,nrci,npc
real(8) cb,t1
! automatic arrays
real(8) rfmt(npcmtmax)
! set the Q=0 muffin-tin potential equal to that of the normal ground-state
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  call rfmtftoc(nrc,nrci,vsmt(:,ias),rfmt)
  call rbshtip(nrc,nrci,rfmt)
  vsqmt(1:npc,ias,1)=rfmt(1:npc)
end do
! zero the muffin-tin potential for non-zero Q
do ifq=2,nfqrz
  do ias=1,natmtot
    is=idxis(ias)
    vsqmt(1:npcmt(is),ias,ifq)=0.d0
  end do
end do
! repeat for the interstitial potential
vsqir(1:ngtot,1)=vsir(1:ngtot)
vsqir(1:ngtot,2:nfqrz)=0.d0
if (.not.spinpol) return
! set the Q=0 muffin-tin magnetic field equal to that of the normal ground-state
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
    bsqmt(1:npc,ias,idm,1)=bsmt(1:npc,ias,idm)
  end do
end do
! zero the magnetic field for non-zero Q
do ifq=2,nfqrz
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      bsqmt(1:npcmt(is),ias,idm,ifq)=0.d0
    end do
  end do
end do
! repeat for the interstitial magnetic field
bsqir(1:ngtot,1:ndmag,1)=bsir(1:ngtot,1:ndmag)
bsqir(1:ngtot,1:ndmag,2:nfqrz)=0.d0
end subroutine

