
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine geomopt
use modmain
use modmpi
use moddelf
use modvars
implicit none
! local variables
integer istp,jstp,is,i
real(8) ds
! initialise global variables (and the muffin-tin radii)
call init0
call init1
! trim the Fourier components of the Kohn-Sham potential
trimvg0=trimvg
trimvg=.true.
! store orginal volume
omega0=omega
! atomic forces are required
tforce=.true.
if (task == 3) then
  trdstate=.true.
else
  trdstate=.false.
end if
! initial atomic step sizes
if (allocated(tauatp)) deallocate(tauatp)
allocate(tauatp(natmtot))
tauatp(:)=tau0atp
! initialise the previous total force on each atom
if (allocated(forcetotp)) deallocate(forcetotp)
allocate(forcetotp(3,natmtot))
forcetotp(:,:)=0.d0
! initial lattice optimisation step size
taulatv(:)=tau0latv
! initialise previous stress tensor
stressp(:)=0.d0
! open formatted files for writing
if (mp_mpi) then
  open(71,file='TOTENERGY_OPT.OUT',form='FORMATTED')
  open(72,file='FORCEMAX.OUT',form='FORMATTED')
  open(73,file='GEOMETRY_OPT.OUT',form='FORMATTED')
  open(74,file='IADIST_OPT.OUT',form='FORMATTED')
  open(75,file='FORCES_OPT.OUT',form='FORMATTED')
  if (spinpol) then
    open(78,file='MOMENTM_OPT.OUT',form='FORMATTED')
  end if
  if (latvopt /= 0) then
    open(86,file='STRESSMAX.OUT',form='FORMATTED')
    open(87,file='STRESS_OPT.OUT',form='FORMATTED')
    open(88,file='OMEGA_OPT.OUT',form='FORMATTED')
  end if
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
do istp=1,maxlatvstp
  do jstp=1,maxatpstp
    if (atpopt == 0) exit
    if (mp_mpi) then
      write(*,'("Info(geomopt): atomic position optimisation step : ",I6)') jstp
    end if
! ground-state and forces calculation
    call gndstate
! check for stop signal
    if (tstop) goto 10
! subsequent calculations will read in the potential from STATE.OUT
    trdstate=.true.
! update the atomic positions
    call atpstep
! write total energy, forces, atomic positions, interatomic distances to file
    if (mp_mpi) then
      write(71,'(G22.12)') engytot
      flush(71)
      write(72,'(G18.10)') forcemax
      flush(72)
      write(73,*); write(73,*)
      write(73,'("! Atomic position optimisation step : ",I6)') jstp
      call writegeom(73)
      flush(73)
      write(74,*); write(74,*)
      write(74,'("Atomic position optimisation step : ",I6)') jstp
      call writeiad(74)
      flush(74)
      write(75,*); write(75,*)
      write(75,'("Atomic position optimisation step : ",I6)') jstp
      call writeforces(75)
      write(75,*)
      write(75,'("Maximum force magnitude over all atoms (target) : ",G18.10,&
       &" (",G18.10,")")') forcemax,epsforce
      flush(75)
      if (spinpol) then
        write(78,'(G22.12)') momtotm
        flush(78)
      end if
    end if
! broadcast forcemax from master process to all other processes
    call mpi_bcast(forcemax,1,mpi_double_precision,0,mpicom,ierror)
! check force convergence
    if (forcemax <= epsforce) then
      if (mp_mpi) then
        write(75,*)
        write(75,'("Force convergence target achieved")')
      end if
      exit
    end if
    if (mp_mpi.and.(jstp == maxatpstp)) then
      write(*,*)
      write(*,'("Warning(geomopt): atomic position optimisation failed to &
       &converge in ",I6," steps")') maxatpstp
    end if
! store the current forces array
    forcetotp(:,:)=forcetot(:,:)
! end loop over atomic position optimisation
  end do
! exit lattice optimisation loop if required
  if (latvopt == 0) exit
  if (mp_mpi) then
    write(*,'("Info(geomopt): lattice optimisation step : ",I6)') istp
  end if
! generate the stress tensor
  call genstress
! take average of current and previous stress tensors
  stress(:)=0.5d0*(stress(:)+stressp(:))
! check for stop signal
  if (tstop) goto 10
! update the lattice vectors
  call latvstep
! write stress tensor components and maximum magnitude to file
  if (mp_mpi) then
    write(71,'(G22.12)') engytot
    flush(71)
    write(73,*); write(73,*)
    write(73,'("! Lattice optimisation step : ",I6)') istp
    call writegeom(73)
    flush(73)
    write(74,*); write(74,*)
    write(74,'("Lattice optimisation step : ",I6)') istp
    call writeiad(74)
    flush(74)
    if (spinpol) then
      write(78,'(G22.12)') momtotm
      flush(78)
    end if
    write(86,'(G18.10)') stressmax
    flush(86)
    write(87,*)
    write(87,'("Lattice optimisation step : ",I6)') istp
    write(87,'("Derivative of total energy w.r.t. strain tensors :")')
    do i=1,nstrain
      write(87,'(G18.10)') stress(i)
    end do
    flush(87)
    write(88,'(G18.10)') omega
    flush(88)
  end if
! check for stress convergence
  if (latvopt == 1) then
    ds=sum(abs(stress(:)))
  else
! stress may be non-zero because of volume constraint; check change in stress
! tensor instead
    ds=sum(abs(stress(:)-stressp(:)))
  end if
! broadcast ds from master process to all other processes
  call mpi_bcast(ds,1,mpi_double_precision,0,mpicom,ierror)
  if ((istp >= 3).and.(ds <= epsstress*tau0latv)) then
    if (mp_mpi) then
      write(87,*)
      write(87,'("Stress convergence target achieved")')
    end if
    exit
  end if
  if (mp_mpi.and.(istp == maxlatvstp)) then
    write(*,*)
    write(*,'("Warning(geomopt): lattice optimisation failed to converge in ",&
     &I6," steps")') maxlatvstp
  end if
  stressp(1:nstrain)=stress(1:nstrain)
! delete the eigenvector files
  call delfiles(evec=.true.)
! end loop over lattice optimisation
end do
10 continue
if (mp_mpi) then
  close(71); close(72); close(73); close(74); close(75)
  if (spinpol) close(78)
  if (latvopt /= 0) then
    close(86); close(87); close(88)
  end if
end if
! ground-state should be run again after lattice optimisation
if (latvopt /= 0) call gndstate
! write optimised variables
if (wrtvars) then
  call writevars('avec (geomopt)',nv=9,rva=avec)
  call writevars('omega (geomopt)',rv=omega)
  do is=1,nspecies
    call writevars('atposl (geomopt)',n1=is,nv=3*natoms(is),rva=atposl(:,:,is))
  end do
  do is=1,nspecies
    call writevars('atposc (geomopt)',n1=is,nv=3*natoms(is),rva=atposc(:,:,is))
  end do
  call writevars('engytot (geomopt)',rv=engytot)
end if
! restore original parameters
trimvg=trimvg0
end subroutine

