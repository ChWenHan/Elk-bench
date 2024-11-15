
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vclcore(wfmt,vmat)
use modmain
use modomp
implicit none
! arguments
complex(4), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: vmat(nstsv,nstsv)
! local variables
integer ist1,ist2,ist3
integer is,ia,ias,m,nthd
integer nrc,nrci,npc
! automatic arrays
complex(4) wfcr(npcmtmax,2),cfmt(npcmtmax)
complex(8) v(nstsv,nstsv)
! allocatable arrays
complex(4), allocatable :: crhomt(:,:)
! external functions
complex(8), external :: zcfmtinp
allocate(crhomt(npcmtmax,nstsv))
call holdthd(nstsv,nthd)
v(1:nstsv,1:nstsv)=0.d0
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist3=1,nstsp(is)
      if (spcore(ist3,is)) then
        do m=-ksp(ist3,is),ksp(ist3,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,ist3,m,npcmtmax,wfcr)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(cfmt,ist1,ist2) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
          do ist1=1,nstsv
! calculate the complex overlap density in spherical harmonics
            if (spinpol) then
              call crho2(npc,wfcr,wfcr(:,2),wfmt(:,ias,1,ist1), &
               wfmt(:,ias,2,ist1),cfmt)
            else
              call crho1(npc,wfcr,wfmt(:,ias,1,ist1),cfmt)
            end if
            call cfsht(nrc,nrci,cfmt,crhomt(:,ist1))
          end do
!$OMP END DO
!$OMP DO
          do ist2=1,nstsv
            call cpotclmt(nrc,nrci,nrcmtmax,rlcmt(:,:,is),wprcmt(:,:,is), &
             crhomt(:,ist2),cfmt)
            do ist1=1,ist2
              v(ist1,ist2)=v(ist1,ist2)-zcfmtinp(nrc,nrci,wrcmt(:,is), &
               crhomt(:,ist1),cfmt)
            end do
          end do
!$OMP END DO
!$OMP END PARALLEL
        end do
      end if
    end do
  end do
end do
call freethd(nthd)
do ist1=1,nstsv
! set the lower triangular part of the matrix
  do ist2=1,ist1-1
    v(ist1,ist2)=conjg(v(ist2,ist1))
  end do
! make the diagonal elements real
  v(ist1,ist1)=dble(v(ist1,ist1))
end do
! scale the Coulomb matrix elements in the case of a hybrid functional
if (hybrid) v(1:nstsv,1:nstsv)=hybridc*v(1:nstsv,1:nstsv)
! add to input matrix
vmat(1:nstsv,1:nstsv)=vmat(1:nstsv,1:nstsv)+v(1:nstsv,1:nstsv)
deallocate(crhomt)
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

