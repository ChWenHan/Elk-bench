
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine moldyn
use modmain
use modtddft
use modmpi
use modomp
use moddelf
implicit none
! local variables
integer is,ia
! initialise universal variables
call init0
! store original parameters
atposl0(:,:,:)=atposl(:,:,:)
atposc0(:,:,:)=atposc(:,:,:)
tshift0=tshift
tforce0=tforce
tfav00=tfav0
trimvg0=trimvg
! no shifting of atomic basis allowed
tshift=.false.
! calculate atomic forces
tforce=.true.
! average force can be non-zero (allow for translation of atomic basis)
tfav0=.false.
! trim the Fourier components of the exchange-correlation potential
trimvg=.true.
! generate the time step grid
call gentimes
if (task == 420) then
! calculation starts at t=0
  tdt0=.true.
else
! calculation restarts
  tdt0=.false.
end if
if (tdt0) then
! start from t=0
  itimes0=1
else
! restart if required
  call readtimes
  itimes0=itimes0+ntsforce
  trdatdv=.true.
end if
if (trdatdv) then
! read the atomic displacements and velocities
  call readatdvc
else
! set the displacements and velocities to zero
  atdvc(:,:,:,:)=0.d0
end if
trdstate=.false.
if (tdt0.and.mp_mpi) then
  call delfile('TOTENERGY_TD.OUT')
  if (spinpol) then
    call delfile('MOMENT_TD.OUT')
    call delfile('MOMENTM_TD.OUT')
    call delfile('MOMENTMT_TD.OUT')
    call delfile('MOMENTIR_TD.OUT')
  end if
  call delfile('FORCETOT_TD.OUT')
  call delfile('FORCEMAX_TD.OUT')
  call delfile('ATDISPL_TD.OUT')
  call delfile('ATDISPC_TD.OUT')
end if
do itimes=itimes0,ntimes-1,ntsforce
  if (mp_mpi) then
    write(*,'("Info(moldyn): time step ",I8," of ",I8,",   t = ",G18.10)') &
     itimes,ntimes,times(itimes)
  end if
! reset the OpenMP thread variables
  call omp_reset
! add the displacements to the atomic positions
  do is=1,nspecies
    do ia=1,natoms(is)
      atposc(:,ia,is)=atposc0(:,ia,is)+atdvc(:,0,ia,is)
      call r3mv(ainv,atposc(:,ia,is),atposl(:,ia,is))
    end do
  end do
! calculate the ground-state and atomic forces
  call gndstate
! subsequent calculations will read in the potential from STATE.OUT
  trdstate=.true.
! time step the atomic positions within the adiabatic approximation
  call atptstep(forcetot)
  if (mp_mpi) then
! write the time step to file
    call writetimes
! write time-dependent total energy
    call writetdengy
! write spin moments if required
    if (spinpol) call writemomtd
! write time-dependent atomic forces
    call writetdforces
! write the time-dependent atomic displacements
    call writeatdisp
! write the XCrysden animation file crystal.axsf
    call writeaxsf
  end if
  if (tstop) exit
end do
! restore original input parameters
atposl(:,:,:)=atposl0(:,:,:)
tshift=tshift0
tforce=tforce0
tfav0=tfav00
trimvg=trimvg0
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

