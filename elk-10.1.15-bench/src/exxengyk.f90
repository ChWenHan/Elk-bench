
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengyk(ikp)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
! local variables
integer iq,ik,jk,m
integer nst1,nst2,ist,jst
integer is,ia,ias
integer nrc,nrci,npc
integer iv(3),ig
real(8) ex,vc(3)
complex(8) z1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqrmt(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(4), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:),wfcr(:,:)
complex(4), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(4), allocatable :: crhomt(:,:),crhoir(:)
complex(4), allocatable :: cvclmt(:,:),cvclir(:)
! external functions
complex(8), external :: zcfinp,zcfmtinp
! get the eigenvectors from file for input reduced k-point
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
call getevecfv(filext,ikp,vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
call getevecsv(filext,ikp,vkl(:,ikp),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk(1,ikp),vgkc(:,:,1,ikp),gkc(:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! count and index the occupied states
nst1=0
do ist=1,nstsv
  if (evalsv(ist,ikp) > efermi) cycle
  nst1=nst1+1
  idx(nst1)=ist
end do
! calculate the wavefunctions for occupied states of the input k-point
allocate(wfmt1(npcmtmax,natmtot,nspinor,nst1),wfir1(ngtc,nspinor,nst1))
call genwfsv_sp(.false.,.false.,nst1,idx,ngdgc,igfc,ngk(1,ikp),igkig(:,1,ikp), &
 apwalm,evecfv,evecsv,wfmt1,ngtc,wfir1)
! allocate local arrays
allocate(vgqc(3,ngvc),gqc(ngvc),gclgq(ngvc))
allocate(jlgqrmt(0:lnpsd,ngvc,nspecies))
allocate(ylmgq(lmmaxo,ngvc),sfacgq(ngvc,natmtot))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv))
allocate(wfir2(ngtc,nspinor,nstsv))
allocate(crhomt(npcmtmax,natmtot),crhoir(ngtc))
allocate(cvclmt(npcmtmax,natmtot),cvclir(ngtc))
! zero the local exchange energy variable
ex=0.d0
! start loop over non-reduced k-point set
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! determine the q-vector
  iv(:)=ivk(:,ikp)-ivk(:,ik)
  iv(:)=modulo(iv(:),ngridk(:))
! check if the q-point is in user-defined set
  iv(:)=iv(:)*ngridq(:)
  if (any(mod(iv(:),ngridk(:)) /= 0)) cycle
  iv(:)=iv(:)/ngridk(:)
  iq=ivqiq(iv(1),iv(2),iv(3))
  vc(:)=vkc(:,ikp)-vkc(:,ik)
  do ig=1,ngvc
! determine the G+q-vectors
    vgqc(1:3,ig)=vgc(1:3,ig)+vc(1:3)
! G+q-vector length
    gqc(ig)=sqrt(vgqc(1,ig)**2+vgqc(2,ig)**2+vgqc(3,ig)**2)
! spherical harmonics for G+q-vectors
    call genylmv(.true.,lmaxo,vgqc(:,ig),ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvc,vgqc,ngvc,sfacgq)
! generate the regularised Coulomb Green's function in G+q-space
  call gengclgq(.true.,iq,ngvc,gqc,gclgq)
! compute the required spherical Bessel functions
  call genjlgprmt(lnpsd,ngvc,gqc,ngvc,jlgqrmt)
! find the matching coefficients
  call match(ngk(1,ik),vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
! get the eigenvectors from file for non-reduced k-point
  call getevecfv(filext,0,vkl(:,ik),vgkl(:,:,1,ik),evecfv)
  call getevecsv(filext,0,vkl(:,ik),evecsv)
! count and index the occupied states
  nst2=0
  do jst=1,nstsv
    if (evalsv(jst,jk) > efermi) cycle
    nst2=nst2+1
    idx(nst2)=jst
  end do
! calculate the wavefunctions for occupied states
  call genwfsv_sp(.false.,.false.,nst2,idx,ngdgc,igfc,ngk(1,ik),igkig(:,1,ik), &
   apwalm,evecfv,evecsv,wfmt2,ngtc,wfir2)
!--------------------------------------------!
!    valence-valence-valence contribution    !
!--------------------------------------------!
  do jst=1,nst2
    do ist=1,nst1
! calculate the complex overlap density
      call gencrho(.true.,.true.,ngtc,wfmt2(:,:,:,jst),wfir2(:,:,jst), &
       wfmt1(:,:,:,ist),wfir1(:,:,ist),crhomt,crhoir)
! calculate the Coulomb potential
      call gencvclmt(nrcmt,nrcmti,nrcmtmax,rlcmt,wprcmt,npcmtmax,crhomt,cvclmt)
      call cpotcoul(nrcmt,nrcmti,npcmt,nrcmtmax,rlcmt,ngdgc,igfc,ngvc,gqc, &
       gclgq,ngvc,jlgqrmt,ylmgq,sfacgq,crhoir,npcmtmax,cvclmt,cvclir)
      cvclir(1:ngtc)=cvclir(1:ngtc)*cfrc(1:ngtc)
      z1=zcfinp(crhomt,crhoir,cvclmt,cvclir)
      ex=ex-0.5d0*occmax*wkpt(ikp)*wqptnr*dble(z1)
    end do
  end do
! end loop over non-reduced k-point set
end do
deallocate(vgqc,gqc,gclgq,jlgqrmt)
deallocate(evecfv,evecsv)
deallocate(apwalm,ylmgq,sfacgq)
deallocate(wfmt2,wfir2)
!-----------------------------------------!
!    valence-core-valence contribution    !
!-----------------------------------------!
allocate(wfcr(npcmtmax,2))
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do jst=1,nstsp(is)
      if (spcore(jst,is)) then
        do m=-ksp(jst,is),ksp(jst,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,jst,m,npcmtmax,wfcr)
          do ist=1,nst1
! calculate the complex overlap density in spherical harmonics
            if (spinpol) then
              call crho2(npc,wfcr,wfcr(:,2),wfmt1(:,ias,1,ist), &
               wfmt1(:,ias,2,ist),crhomt(:,ias))
            else
              call crho1(npc,wfcr,wfmt1(:,ias,1,ist),crhomt(:,ias))
            end if
            call cfshtip(nrc,nrci,crhomt(:,ias))
! calculate the Coulomb potential
            call cpotclmt(nrc,nrci,nrcmtmax,rlcmt(:,:,is),wprcmt(:,:,is), &
             crhomt(:,ias),cvclmt(:,ias))
            z1=zcfmtinp(nrc,nrci,wrcmt(:,is),crhomt(:,ias),cvclmt(:,ias))
            ex=ex-occmax*wkpt(ikp)*dble(z1)
          end do
! end loop over m
        end do
! end loop over jst
      end if
    end do
! end loops over atoms and species
  end do
end do
! add to global exchange energy
!$OMP CRITICAL(exxengyk_)
engyx=engyx+ex
!$OMP END CRITICAL(exxengyk_)
deallocate(wfmt1,wfir1,wfcr)
deallocate(crhomt,crhoir,cvclmt,cvclir)
return

contains

pure subroutine crho1(n,wf1,wf2,crho)
implicit none
integer, intent(in) :: n
complex(4), intent(in) :: wf1(n),wf2(n)
complex(4), intent(out) :: crho(n)
crho(1:n)=conjg(wf1(1:n))*wf2(1:n)
end subroutine

pure subroutine crho2(n,wf11,wf12,wf21,wf22,crho)
implicit none
integer, intent(in) :: n
complex(4), intent(in) :: wf11(n),wf12(n),wf21(n),wf22(n)
complex(4), intent(out) :: crho(n)
crho(1:n)=conjg(wf11(1:n))*wf21(1:n)+conjg(wf12(1:n))*wf22(1:n)
end subroutine

end subroutine

