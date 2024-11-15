
! Copyright (C) 2015 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writetddos
! !INTERFACE:
subroutine writetddos
! !USES:
use modmain
use modtddft
use modmpi
! !DESCRIPTION:
!   Calculates the time-dependent density of states (DOS). This is defined as
!   $$ {\rm DOS}(\omega,t)=\frac{\Omega}{(2\pi)^3}\int d^3k \sum_i
!    \delta(\varepsilon_{i{\bf k}}-\omega) F_{i{\bf k}}(t), $$
!   where
!   $$ F_{i{\bf k}}(t)=\sum_j f_{j{\bf k}}|\langle\varphi_{i{\bf k}}|
!    \phi_{j{\bf k}}(t)\rangle|^2, $$
!   with occupation numbers $f_{j{\bf k}}$, ground-state orbitals
!   $\varphi_{i{\bf k}}$ and time-dependent orbitals $\phi_{j{\bf k}}(t)$.
!
! !REVISION HISTORY:
!   Created April 2015 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,jst,lp
real(8) sm,t1
complex(8) z1
character(256) fext
! allocatable arrays
real(8), allocatable :: occsvp(:,:)
complex(8), allocatable :: evecsv(:,:),evecsvt(:,:)
! external functions
complex(8), external :: zdotc
! file extension
write(fext,'("_TS",I8.8,".OUT")') itimes
allocate(occsvp(nstsv,nkpt))
allocate(evecsv(nstsv,nstsv),evecsvt(nstsv,nstsv))
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
! read in ground-state eigenvectors
  call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! read in the time evolving eigenvectors
  call getevecsv('_TD.OUT',ik,vkl(:,ik),evecsvt)
! determine the time-dependent projected occupation numbers
  do ist=1,nstsv
    sm=0.d0
    do jst=1,nstsv
      t1=occsv(jst,ik)
      if (abs(t1) < epsocc) cycle
      z1=zdotc(nstsv,evecsv(:,ist),1,evecsvt(:,jst),1)
      sm=sm+t1*(dble(z1)**2+aimag(z1)**2)
    end do
    occsvp(ist,ik)=sm
  end do
! write projected occupation numbers to file
  call putoccsv('P'//trim(fext),ik,occsvp(:,ik))
end do
deallocate(evecsv,evecsvt)
! broadcast projected occupation numbers to every MPI process
if (np_mpi > 1) then
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(occsvp(:,ik),nstsv,mpi_double_precision,lp,mpicom,ierror)
  end do
end if
if (mp_mpi) then
! compute the effective electronic temperature and write to file
  call tdtemp(occsvp)
! write the DOS to file
  call dos(fext,.true.,occsvp)
end if
deallocate(occsvp)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine
!EOC

