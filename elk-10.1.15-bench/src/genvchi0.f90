
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvchi0(t3hw,ik,lock,vqpl,gclgq,jlgqr,ylmgq,sfacgq,nm,vchi0)
use modmain
use modomp
implicit none
! local variables
logical, intent(in) :: t3hw
integer, intent(in) :: ik
integer(omp_lock_kind), intent(inout) :: lock(nwrf)
real(8), intent(in) :: vqpl(3),gclgq(ngrf)
real(8), intent(in) :: jlgqr(njcmax,nspecies,ngrf)
complex(8), intent(in) :: ylmgq(lmmaxo,ngrf)
complex(8), intent(in) :: sfacgq(ngrf,natmtot)
integer, intent(in) :: nm
complex(8), intent(inout) :: vchi0(nm,nm,nwrf)
! local variables
logical tq0
integer isym,jk,jkq,iw
integer nst,nstq,ist,jst,kst,lst
integer nm2,ig,jg,i,j,nthd
real(8) vkql(3),ei,ej,eij,t1
complex(8) a(3,3),z1
! automatic arrays
integer idx(nstsv),idxq(nstsv)
integer ngp(nspnfv),ngpq(nspnfv)
! allocatable arrays
integer, allocatable :: igpig(:,:),igpqig(:,:)
complex(4), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(4), allocatable :: wfmtq(:,:,:,:),wfirq(:,:,:)
complex(4), allocatable :: crhomt(:,:),crhoir(:)
complex(8), allocatable :: zrhoig(:),pmat(:,:,:),zw(:),b(:,:)
! check if q = 0
tq0=.false.
if (sum(abs(vqpl(:))) < epslat) tq0=.true.
! k+q-vector in lattice coordinates
vkql(:)=vkl(:,ik)+vqpl(:)
! equivalent reduced k-points for k and k+q
jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
call findkpt(vkql,isym,jkq)
! count and index states at k and k+q in energy window
nst=0
do ist=1,nstsv
  if (abs(evalsv(ist,jk)-efermi) > emaxrf) cycle
  nst=nst+1
  idx(nst)=ist
end do
nstq=0
do ist=1,nstsv
  if (abs(evalsv(ist,jkq)-efermi) > emaxrf) cycle
  nstq=nstq+1
  idxq(nstq)=ist
end do
! generate the wavefunctions for all states at k and k+q in energy window
allocate(igpig(ngkmax,nspnfv))
allocate(wfmt(npcmtmax,natmtot,nspinor,nst),wfir(ngtc,nspinor,nst))
call genwfsvp_sp(.false.,.false.,nst,idx,ngdgc,igfc,vkl(:,ik),ngp,igpig,wfmt, &
 ngtc,wfir)
deallocate(igpig)
allocate(igpqig(ngkmax,nspnfv))
allocate(wfmtq(npcmtmax,natmtot,nspinor,nstq),wfirq(ngtc,nspinor,nstq))
call genwfsvp_sp(.false.,.false.,nstq,idxq,ngdgc,igfc,vkql,ngpq,igpqig,wfmtq, &
 ngtc,wfirq)
deallocate(igpqig)
! read the momentum matrix elements from file for q = 0
if (tq0) then
  allocate(pmat(nstsv,nstsv,3))
  call getpmat(vkl(:,ik),pmat)
! divide by unit cell volume
  t1=1.d0/omega
  pmat(1:nstsv,1:nstsv,1:3)=t1*pmat(1:nstsv,1:nstsv,1:3)
end if
nm2=nm**2
call holdthd(nst,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(crhomt,crhoir,zrhoig,zw,b) &
!$OMP PRIVATE(jst,kst,lst,ei,ej,eij,t1) &
!$OMP PRIVATE(iw,ig,jg,z1,i,j,a) &
!$OMP NUM_THREADS(nthd)
allocate(crhomt(npcmtmax,natmtot),crhoir(ngtc))
allocate(zrhoig(ngrf),zw(nwrf))
if (tq0.and.t3hw) then
  allocate(b(-1:ngrf,-1:ngrf))
else
  allocate(b(ngrf,ngrf))
end if
!$OMP DO
do ist=1,nst
  kst=idx(ist)
  ei=evalsv(kst,jk)
  do jst=1,nstq
    lst=idxq(jst)
    t1=wkptnr*omega*(occsv(kst,jk)-occsv(lst,jkq))
    if (abs(t1) < 1.d-8) cycle
    ej=evalsv(lst,jkq)
    eij=ei-ej
! frequency-dependent part in response function formula for all frequencies
    do iw=1,nwrf
      zw(iw)=t1/(eij+wrf(iw))
    end do
! compute the complex density in G+q-space
    call gencrho(.true.,.true.,ngtc,wfmt(:,:,:,ist),wfir(:,:,ist), &
     wfmtq(:,:,:,jst),wfirq(:,:,jst),crhomt,crhoir)
    call zftcf(ngrf,jlgqr,ylmgq,ngrf,sfacgq,crhomt,crhoir,zrhoig)
! Hermitian part of body
    do jg=1,ngrf
      do ig=1,jg-1
        b(ig,jg)=conjg(b(jg,ig))
      end do
      z1=gclgq(jg)*conjg(zrhoig(jg))
      do ig=jg,ngrf
        b(ig,jg)=gclgq(ig)*zrhoig(ig)*z1
      end do
    end do
! case of q = 0
    if (tq0) then
      if (t3hw) then
        b(-1:1,-1:1)=0.d0
! calculate 3 x ngrf wings of matrix
        t1=-sqrt(fourpi)/eij
        do i=-1,1
          z1=t1*pmat(kst,lst,i+2)
          b(i,2:ngrf)=z1*conjg(zrhoig(2:ngrf))*gclgq(2:ngrf)
          do j=2,ngrf
            b(j,i)=conjg(b(i,j))
          end do
        end do
      else
! use trace of 3 x 3 head of matrix
        t1=sum(dble(pmat(kst,lst,1:3))**2+aimag(pmat(kst,lst,1:3))**2)/3.d0
        b(1,1)=(fourpi/eij**2)*t1
! wings of matrix
        t1=-sqrt(fourpi)/eij
        z1=(t1/3.d0)*(pmat(kst,lst,1)+pmat(kst,lst,2)+pmat(kst,lst,3))
        b(1,2:ngrf)=z1*conjg(zrhoig(2:ngrf))*gclgq(2:ngrf)
        b(2:ngrf,1)=conjg(b(1,2:ngrf))
      end if
    end if
! add to body and wings of the response function
    do iw=1,nwrf
      call omp_set_lock(lock(iw))
      call zaxpy(nm2,zw(iw),b,1,vchi0(1,1,iw),1)
      call omp_unset_lock(lock(iw))
    end do
! calculate 3 x 3 head
    if (tq0.and.t3hw) then
      t1=-fourpi/eij
      zw(1:nwrf)=zw(1:nwrf)/wrf(1:nwrf)
      do j=1,3
        do i=1,3
          a(i,j)=t1*pmat(kst,lst,i)*conjg(pmat(kst,lst,j))
        end do
      end do
! add to the head of the response function
      do iw=1,nwrf
        call omp_set_lock(lock(iw))
        vchi0(1:3,1:3,iw)=vchi0(1:3,1:3,iw)+a(1:3,1:3)*zw(iw)
        call omp_unset_lock(lock(iw))
      end do
    end if
! end loop over jst
  end do
! end loop over ist
end do
!$OMP END DO
deallocate(crhomt,crhoir,zrhoig,zw,b)
!$OMP END PARALLEL
call freethd(nthd)
deallocate(wfmt,wfir,wfmtq,wfirq)
if (tq0) deallocate(pmat)
end subroutine

