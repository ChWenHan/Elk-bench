
! Copyright (C) 2006 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnhf(ikp,vmt,vir,bmt,bir,evecsvp)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtc)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtc,ndmag)
complex(8), intent(inout) :: evecsvp(nstsv,nstsv)
! local variables
integer ik,jk,nst
integer ist1,ist2,ist3,jst3
integer iv(3),iq,ig,nthd
real(8) vc(3),t1
complex(8) z1
! automatic arrays
integer idx(nstsv)
real(8) vgqc(3,ngvc),gqc(ngvc),gclgq(ngvc)
! allocatable arrays
real(8), allocatable :: jlgqrmt(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: h(:,:),v(:,:),kmat(:,:)
complex(4), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:)
complex(4), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(4), allocatable :: crhomt(:,:,:),crhoir(:,:)
complex(4), allocatable :: cvclmt(:,:),cvclir(:)
! external functions
complex(8), external :: zcfinp
!$OMP CRITICAL(eveqnhf_)
write(*,'("Info(eveqnhf): ",I6," of ",I6," k-points")') ikp,nkpt
!$OMP END CRITICAL(eveqnhf_)
! allocate local arrays
allocate(jlgqrmt(0:lnpsd,ngvc,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxo,ngvc),sfacgq(ngvc,natmtot))
allocate(h(nstsv,nstsv),v(nstsv,nstsv))
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv),wfir1(ngtc,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv),wfir2(ngtc,nspinor,nstsv))
allocate(crhomt(npcmtmax,natmtot,nstsv),crhoir(ngtc,nstsv))
! get the first-variational eigenvectors from file for input reduced k-point
call getevecfv(filext,ikp,vkl(:,ikp),vgkl(:,:,1,ikp),evecfv)
! find the matching coefficients
call match(ngk(1,ikp),vgkc(:,:,1,ikp),gkc(:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! calculate the wavefunctions for all states of the input k-point
call genwfsv_sp(.false.,.true.,nstsv,[0],ngdgc,igfc,ngk(1,ikp),igkig(:,1,ikp), &
 apwalm,evecfv,evecsvp,wfmt1,ngtc,wfir1)
!-----------------------------------------!
!     local potential matrix elements     !
!-----------------------------------------!
if (hybrid.and.spinpol) then
! magnetic field matrix elements in hybrid case
  call genvbmatk(vmt,vir,bmt,bir,ngk(1,ikp),igkig(:,1,ikp),wfmt1,ngtc,wfir1,h)
else
  call genvmatk(vmt,vir,ngk(1,ikp),igkig(:,1,ikp),wfmt1,ngtc,wfir1,h)
end if
! Fourier transform wavefunctions to real-space
call cftwfir(ngk(1,ikp),igkig(:,1,ikp),wfir1)
!---------------------------------!
!     kinetic matrix elements     !
!---------------------------------!
allocate(kmat(nstsv,nstsv))
call getkmat(ikp,kmat)
call zgemm('N','N',nstsv,nstsv,nstsv,zone,kmat,nstsv,evecsvp,nstsv,zzero,v, &
 nstsv)
call zgemm('C','N',nstsv,nstsv,nstsv,zone,evecsvp,nstsv,v,nstsv,zone,h,nstsv)
deallocate(kmat)
!------------------------------!
!     Fock matrix elements     !
!------------------------------!
v(:,:)=0.d0
! loop over non-reduced k-point set
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
! count and index occupied states
  nst=0
  do ist3=1,nstsv
    if (abs(occsv(ist3,jk)) < epsocc) cycle
    nst=nst+1
    idx(nst)=ist3
  end do
! calculate the wavefunctions for occupied states
  call genwfsv_sp(.false.,.false.,nst,idx,ngdgc,igfc,ngk(1,ik),igkig(:,1,ik), &
   apwalm,evecfv,evecsv,wfmt2,ngtc,wfir2)
  call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(cvclmt,cvclir) &
!$OMP PRIVATE(ist1,ist2,ist3,jst3,t1,z1) &
!$OMP NUM_THREADS(nthd)
  allocate(cvclmt(npcmtmax,natmtot),cvclir(ngtc))
  do ist3=1,nst
    jst3=idx(ist3)
! calculate the complex overlap densities for all states (T. McQueen)
!$OMP DO
    do ist1=1,nstsv
      call gencrho(.true.,.true.,ngtc,wfmt2(:,:,:,ist3),wfir2(:,:,ist3), &
       wfmt1(:,:,:,ist1),wfir1(:,:,ist1),crhomt(:,:,ist1),crhoir(:,ist1))
    end do
!$OMP END DO
    t1=wqptnr*occsv(jst3,jk)/occmax
!$OMP DO
    do ist2=1,nstsv
! calculate the Coulomb potential
      call gencvclmt(nrcmt,nrcmti,nrcmtmax,rlcmt,wprcmt,npcmtmax, &
       crhomt(:,:,ist2),cvclmt)
      call cpotcoul(nrcmt,nrcmti,npcmt,nrcmtmax,rlcmt,ngdgc,igfc,ngvc,gqc, &
       gclgq,ngvc,jlgqrmt,ylmgq,sfacgq,crhoir(:,ist2),npcmtmax,cvclmt,cvclir)
      cvclir(:)=cvclir(:)*cfrc(:)
      do ist1=1,ist2
        z1=zcfinp(crhomt(:,:,ist1),crhoir(:,ist1),cvclmt,cvclir)
        v(ist1,ist2)=v(ist1,ist2)-t1*z1
      end do
    end do
!$OMP END DO
  end do
  deallocate(cvclmt,cvclir)
!$OMP END PARALLEL
  call freethd(nthd)
! end loop over non-reduced k-point set
end do
deallocate(jlgqrmt,ylmgq,sfacgq,apwalm,evecfv)
deallocate(wfmt1,wfir1,wfmt2,wfir2,crhomt,crhoir)
! scale the Coulomb matrix elements in the case of a hybrid functional
if (hybrid) v(:,:)=hybridc*v(:,:)
! add the Coulomb matrix elements to Hamiltonian
h(:,:)=h(:,:)+v(:,:)
!----------------------------------------------!
!     diagonalise Hartree-Fock Hamiltonian     !
!----------------------------------------------!
call eveqnzh(nstsv,nstsv,h,evalsv(:,ikp))
! apply unitary transformation to the third-variational states so that they
! refer to the first-variational basis
evecsv(:,:)=evecsvp(:,:)
call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,h,nstsv,zzero,evecsvp, &
 nstsv)
deallocate(evecsv,h,v)
end subroutine

