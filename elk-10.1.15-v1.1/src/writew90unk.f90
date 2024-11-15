
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writew90unk
use modmain
use modw90
use modomp
implicit none
! local variables
integer ispn,ik,ist,nthd
integer is,ias,nrc,nrci,npc
integer np,ngp,iu,i
real(8) vc(3),kc
character(256) fname
! automatic arrays
complex(8) ylmk(lmmaxo),sfack(natmtot)
! allocatable arrays
integer, allocatable :: igpig(:,:)
real(8), allocatable :: vpl(:,:),jlkr(:,:)
real(8), allocatable :: wfmta(:,:),wfmtb(:,:)
real(8), allocatable :: wfira(:),wfirb(:)
real(8), allocatable :: wfa(:,:,:),wfb(:,:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:),expmt(:,:)
! total number of plot points
np=np3d(1)*np3d(2)*np3d(3)
! generate the 3D plotting points
allocate(vpl(3,np))
call plotpt3d(vpl)
! parallel loop over non-reduced k-points
call holdthd(nkptnr,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(igpig,jlkr,wfmta,wfmtb) &
!$OMP PRIVATE(wfira,wfirb,wfa,wfb) &
!$OMP PRIVATE(wfmt,wfir,ylmk,sfack,expmt) &
!$OMP PRIVATE(ngp,vc,kc,ist,ispn,ias,is) &
!$OMP PRIVATE(nrc,nrci,npc,fname,iu,i) &
!$OMP NUM_THREADS(nthd)
allocate(igpig(ngkmax,nspnfv))
allocate(jlkr(njcmax,nspecies))
allocate(wfmta(npmtmax,natmtot),wfmtb(npmtmax,natmtot))
allocate(wfira(ngtot),wfirb(ngtot))
allocate(wfa(np,nspinor,num_bands),wfb(np,nspinor,num_bands))
allocate(wfmt(npcmtmax,natmtot,nspinor,num_bands))
allocate(wfir(ngtot,nspinor,num_bands))
allocate(expmt(npcmtmax,natmtot))
!$OMP DO
do ik=1,nkptnr
!$OMP CRITICAL(writew90unk_)
  write(*,'("Info(writew90unk): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL(writew90unk_)
! generate the second-variational wavefunctions
  call genwfsvp(.false.,.false.,num_bands,idxw90,ngridg,igfft,vkl(:,ik),ngp, &
   igpig,wfmt,ngtot,wfir)
! generate the phase factor function exp(-ik.r) in the muffin-tins
  vc(:)=-vkc(:,ik)
  kc=sqrt(vc(1)**2+vc(2)**2+vc(3)**2)
  call genjlgpr(1,kc,jlkr)
  call genylmv(.true.,lmaxo,vc,ylmk)
  call gensfacgp(1,vc,1,sfack)
  call genexpmt(1,jlkr,ylmk,1,sfack,expmt)
! split the wavefunctions into real and imaginary parts
  do ist=1,num_bands
    do ispn=1,nspinor
      do ias=1,natmtot
        is=idxis(ias)
        nrc=nrcmt(is)
        nrci=nrcmti(is)
        npc=npcmt(is)
! remove the explicit phase exp(ik.r) from the muffin-tin wavefunction
        wfmt(1:npc,ias,ispn,ist)=wfmt(1:npc,ias,ispn,ist)*expmt(1:npc,ias)
        wfmta(1:npc,ias)=dble(wfmt(1:npc,ias,ispn,ist))
        call rfshtip(nrc,nrci,wfmta(:,ias))
        wfmtb(1:npc,ias)=aimag(wfmt(1:npc,ias,ispn,ist))
        call rfshtip(nrc,nrci,wfmtb(:,ias))
      end do
      call rfmtctof(wfmta)
      call rfmtctof(wfmtb)
      wfira(:)=dble(wfir(:,ispn,ist))
      wfirb(:)=aimag(wfir(:,ispn,ist))
! generate the wavefunctions on a regular grid
      call rfplot(np,vpl,wfmta,wfira,wfa(:,ispn,ist))
      call rfplot(np,vpl,wfmtb,wfirb,wfb(:,ispn,ist))
    end do
  end do
  if (spinpol) then
    write(fname,'("UNK",I5.5,".NC")') ik
  else
    write(fname,'("UNK",I5.5,".1")') ik
  end if
  open(newunit=iu,file=trim(fname),form='UNFORMATTED',action='WRITE')
  write(iu) np3d(1),np3d(2),np3d(3),ik,num_bands
  do ist=1,num_bands
    write(iu) (cmplx(wfa(i,1,ist),wfb(i,1,ist),8),i=1,np)
    if (spinpol) then
      write(iu) (cmplx(wfa(i,2,ist),wfb(i,2,ist),8),i=1,np)
    end if
  end do
  close(iu)
end do
!$OMP END DO
deallocate(igpig,jlkr)
deallocate(wfmta,wfmtb,wfira,wfirb)
deallocate(wfa,wfb,wfmt,wfir,expmt)
!$OMP END PARALLEL
call freethd(nthd)
write(*,*)
write(*,'("Info(writew90unk): created the UNKkkkkk.s files")')
deallocate(vpl)
end subroutine

