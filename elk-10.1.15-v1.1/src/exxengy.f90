
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengy
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,ist,jst,is,ia
integer nrc,nrci,npc
integer m1,m2,nthd
complex(8) z1
! automatic arrays
complex(4) wfcr1(npcmtmax,2),wfcr2(npcmtmax,2)
complex(4) crhomt(npcmtmax),cvclmt(npcmtmax)
! external functions
complex(8), external :: zcfmtinp
! zero the exchange energy
engyx=0.d0
!--------------------------------------------------!
!     val-val-val and val-cr-val contributions     !
!--------------------------------------------------!
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
!$OMP CRITICAL(exxengy_)
  write(*,'("Info(exxengy): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(exxengy_)
  call exxengyk(ik)
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
! add energies from each process and redistribute
call mpi_allreduce(mpi_in_place,engyx,1,mpi_double_precision,mpi_sum,mpicom, &
 ierror)
!-----------------------------------!
!    core-core-core contribution    !
!-----------------------------------!
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    do jst=1,nstsp(is)
      if (spcore(jst,is)) then
        do m2=-ksp(jst,is),ksp(jst,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,jst,m2,npcmtmax,wfcr2)
          do ist=1,nstsp(is)
            if (spcore(ist,is)) then
              do m1=-ksp(ist,is),ksp(ist,is)-1
                call wavefcr(.false.,lradstp,is,ia,ist,m1,npcmtmax,wfcr1)
! calculate the complex overlap density
                call crho2(npc,wfcr1,wfcr1(:,2),wfcr2,wfcr2(:,2),crhomt)
                call cfshtip(nrc,nrci,crhomt)
! calculate the Coulomb potential
                call cpotclmt(nrc,nrci,nrcmtmax,rlcmt(:,:,is),wprcmt(:,:,is), &
                 crhomt,cvclmt)
                z1=zcfmtinp(nrc,nrci,wrcmt(:,is),crhomt,cvclmt)
                engyx=engyx-0.5d0*dble(z1)
              end do
! end loop over ist
            end if
          end do
        end do
! end loop over jst
      end if
    end do
! end loops over atoms and species
  end do
end do
return

contains

pure subroutine crho2(n,wf11,wf12,wf21,wf22,crho)
implicit none
integer, intent(in) :: n
complex(4), intent(in) :: wf11(n),wf12(n),wf21(n),wf22(n)
complex(4), intent(out) :: crho(n)
crho(1:n)=conjg(wf11(1:n))*wf21(1:n)+conjg(wf12(1:n))*wf22(1:n)
end subroutine

end subroutine

