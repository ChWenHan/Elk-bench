
! Copyright (C) 2018 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gwrhomag
use modmain
use modgw
use modmpi
use modomp
implicit none
! local variables
integer ik,lp,nthd
! allocatable arrays
real(8), allocatable :: vmt(:,:),vir(:)
real(8), allocatable :: bmt(:,:,:),bir(:,:)
complex(8), allocatable :: se(:,:,:)
! generate the momentum matrix elements
call genpmat
! generate the inverse RPA response function
call epsinv
! compute the matrix elements of -V_xc and -B_xc
allocate(vmt(npcmtmax,natmtot),vir(ngtc))
if (spinpol) then
  allocate(bmt(npcmtmax,natmtot,ndmag),bir(ngtc,ndmag))
end if
call gwlocal(vmt,vir,bmt,bir)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
if (mp_mpi) write(*,*)
! loop over reduced k-point set
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(se) &
!$OMP NUM_THREADS(nthd)
allocate(se(nstsv,nstsv,0:nwfm))
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
!$OMP CRITICAL(gwrhomag_)
  write(*,'("Info(gwrhomag): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(gwrhomag_)
! determine the self-energy at the fermionic frequencies for current k-point
  call gwsefmk(ik,vmt,vir,bmt,bir,se)
! write the self-energy to file
  call putgwsefm(ik,se)
end do
!$OMP END DO
deallocate(se)
!$OMP END PARALLEL
call freethd(nthd)
deallocate(vmt,vir)
if (spinpol) deallocate(bmt,bir)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! determine the GW Fermi energy
call gwefermi
! compute the GW density matrices and write the natural orbitals and occupation
! numbers to EVECSV.OUT and OCCSV.OUT, respectively
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
  call gwdmatk(ik)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! broadcast occupation number array to every MPI process
if (np_mpi > 1) then
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(occsv(:,ik),nstsv,mpi_double_precision,lp,mpicom,ierror)
  end do
end if
! write the occupation numbers to file
if (mp_mpi) then
  do ik=1,nkpt
    call putoccsv(filext,ik,occsv(:,ik))
  end do
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! determine the density and magnetisation in the usual way
call rhomag
end subroutine

