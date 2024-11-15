
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hflocal(vmt,vir,bmt,bir)
use modmain
implicit none
! arguments
real(8), intent(out) :: vmt(npcmtmax,natmtot),vir(ngtc)
real(8), intent(out) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtc,ndmag)
! local variables
integer idm,is,ias
integer np,nrc,nrci
! automatic arrays
real(8) rfmt(npmtmax),rfir(ngtot)
! compute the Coulomb potential
call potcoul
! convert to spherical coordinates and store in output arrays
if (hybrid) then
! hybrid functional case
  call potxc(.true.,xctype,rhomt,rhoir,magmt,magir,taumt,tauir,exmt,exir,ecmt, &
   ecir,vxcmt,vxcir,bxcmt,bxcir,wxcmt,wxcir)
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    np=npmt(is)
    rfmt(1:np)=vclmt(1:np,ias)+vxcmt(1:np,ias)
    call rfmtftoc(nrc,nrci,rfmt,vmt(:,ias))
    call rbshtip(nrc,nrci,vmt(:,ias))
    call rfcmtwr(nrc,nrci,wrcmt(:,is),vmt(:,ias))
  end do
  rfir(1:ngtot)=vclir(1:ngtot)+vxcir(1:ngtot)
! multiply interstitial potential by characteristic function and convert to
! coarse grid
  call rfirftoc(rfir,vir)
  if (spinpol) then
    do idm=1,ndmag
      do ias=1,natmtot
        is=idxis(ias)
        nrc=nrcmt(is)
        nrci=nrcmti(is)
        call rfmtftoc(nrc,nrci,bxcmt(:,ias,idm),bmt(:,ias,idm))
        call rbshtip(nrc,nrci,bmt(:,ias,idm))
        call rfcmtwr(nrc,nrci,wrcmt(:,is),bmt(:,ias,idm))
      end do
      call rfirftoc(bxcir(:,idm),bir(:,idm))
    end do
  end if
else
! normal Hartree-Fock case
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    call rfmtftoc(nrc,nrci,vclmt(:,ias),vmt(:,ias))
    call rbshtip(nrc,nrci,vmt(:,ias))
    call rfcmtwr(nrc,nrci,wrcmt(:,is),vmt(:,ias))
  end do
  call rfirftoc(vclir,vir)
end if
end subroutine

