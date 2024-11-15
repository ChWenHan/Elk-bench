
! Copyright (C) 2018 P. Elliott, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gwefermi
use modmain
use modgw
use modmpi
use modomp
implicit none
! local variables
logical done
integer, parameter :: maxit=1000
integer ik,ist,it,nthd
real(8) e0,e1,e,chg,chgk
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(gwefermi): finding the GW Fermi energy")')
end if
! find minimum and maximum eigenvalues
e0=evalsv(1,1)
e1=e0
do ik=1,nkpt
  do ist=1,nstsv
    e=evalsv(ist,ik)
    if (e < e0) e0=e
    if (e > e1) e1=e
  end do
end do
done=.false.
do it=1,maxit
  if (mp_mpi.and.(mod(it,10) == 0)) then
    write(*,'("Info(gwefermi): done ",I4," iterations")') it
  end if
  efermi=0.5d0*(e0+e1)
  chg=0.d0
! begin parallel loop over k-points
  call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(chgk) REDUCTION(+:chg) &
!$OMP NUM_THREADS(nthd)
  do ik=1,nkpt
! distribute among MPI processes
    if (mod(ik-1,np_mpi) /= lp_mpi) cycle
    call gwchgk(ik,chgk)
    chg=chg+chgk
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
! add charge from each process and redistribute
  if (np_mpi > 1) then
    call mpi_allreduce(mpi_in_place,chg,1,mpi_double_precision,mpi_sum,mpicom, &
     ierror)
  end if
  if (chg < chgval) then
    e0=efermi
  else
    e1=efermi
  end if
! check for convergence
  if ((e1-e0) < 1.d-12) done=.true.
! broadcast done from master process to all other processes
  call mpi_bcast(done,1,mpi_logical,0,mpicom,ierror)
  if (done) return
end do
write(*,*)
write(*,'("Warning(gwefermi): could not find GW Fermi energy")')
end subroutine

