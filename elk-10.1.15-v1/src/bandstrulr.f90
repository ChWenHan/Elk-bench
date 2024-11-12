
! Copyright (C) 2024 Wenhan Chen, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bandstrulr
use modmain
use modulr
use modmpi
use modomp
implicit none
! local variables
integer ik0,ist,iv,lp,nthd
real(8) emin,emax
! allocatable arrays
complex(8), allocatable :: evecu(:,:)
! initialise global variables
call init0
call init1
call readstate
call genvsig
call gencore
call linengy
call genapwlofr
call gensocfr
call genevfsv
call occupy
call initulr
! read in the potential STATE_ULR.OUT
call readstulr
! initialise the external Coulomb potential
call vclqinit
! apply required local operations to the potential and magnetic field
call vblocalu
! loop over original k-points
call holdthd(nkpt0/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecu) &
!$OMP NUM_THREADS(nthd)
allocate(evecu(nstulr,nstulr))
!$OMP DO
do ik0=1,nkpt0
! distribute among MPI processes
  if (mod(ik0-1,np_mpi) /= lp_mpi) cycle
!$OMP CRITICAL(bandstrulr_)
  write(*,'("Info(bandstrulr): ",I6," of ",I6," k-points")') ik0,nkpt0
!$OMP END CRITICAL(bandstrulr_)
! solve the ultra long-range eigenvalue equation
  call eveqnulr(ik0,evecu)
end do
!$OMP END DO
deallocate(evecu)
!$OMP END PARALLEL
call freethd(nthd)
! broadcast eigenvalue array to every process
if (np_mpi > 1) then
  do ik0=1,nkpt0
    lp=mod(ik0-1,np_mpi)
    call mpi_bcast(evalu(:,ik0),nstulr,mpi_double_precision,lp,mpicom,ierror)
  end do
end if
! subtract the Fermi energy
evalu(:,:)=evalu(:,:)-efermi
! find the minimum and maximum eigenvalues
emin=minval(evalu(:,:))
emax=maxval(evalu(:,:))
if (mp_mpi) then
! output the band structure
  open(50,file='BANDULR.OUT',form='FORMATTED',action='WRITE')
  do ist=1,nstulr
    do ik0=1,nkpt0
      write(50,'(2G18.10)') dpp1d(ik0),evalu(ist,ik0)
    end do
    write(50,*)
  end do
  close(50)
! output the vertex location lines
  open(50,file='BANDLINES.OUT',form='FORMATTED',action='WRITE')
  do iv=1,nvp1d
    write(50,'(2G18.10)') dvp1d(iv),emin
    write(50,'(2G18.10)') dvp1d(iv),emax
    write(50,*)
  end do
  close(50)
  write(*,*)
  write(*,'("Info(bandstrulr):")')
  write(*,'(" Ultra long-range band structure plot written to BANDULR.OUT")')
  write(*,*)
  write(*,'(" Vertex location lines written to BANDLINES.OUT")')
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

