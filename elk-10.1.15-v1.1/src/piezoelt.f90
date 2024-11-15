
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine piezoelt
use modmain
use modmpi
use modtest
implicit none
! local variables
integer i,j,k
real(8) pvl1(3),pvl2(3)
real(8) vc(3),t1
! allocatable arrays
real(8), allocatable :: pelt(:,:)
! initialise universal variables
call init0
call init1
! store original parameters
avec0(:,:)=avec(:,:)
ngridk0(:)=ngridk(:)
maxscl0=maxscl
tshift0=tshift
! no shifting of the atomic basis
tshift=.false.
! generate the strain tensors
call genstrain
! allocate the piezoelectric tensor array
allocate(pelt(3,nstrain))
! initial ground-state run should start from atomic densities
trdstate=.false.
! run the ground-state calculation
call gndstate
! subsequent calculations will read in the previous potential
trdstate=.true.
! compute the first polarisation in lattice coordinates
call polar(pvl1)
! loop over strain tensors
do istrain=1,nstrain
  if (mp_mpi) then
    write(*,'("Info(piezoelt): working on strain tensor ",I1," of ",I1)') &
     istrain,nstrain
  end if
! run the ground-state calculation again
  call gndstate
! compute the second polarisation
  call polar(pvl2)
  do i=1,3
! add multiple of 2*pi to bring polarisation vectors into coincidence
    pvl1(i)=modulo(pvl1(i),twopi)
    pvl2(i)=modulo(pvl2(i),twopi)
    t1=pvl1(i)-pvl2(i)
    if (abs(t1-twopi) < abs(t1)) then
      pvl1(i)=pvl1(i)-twopi
    else if (abs(t1+twopi) < abs(t1)) then
      pvl1(i)=pvl1(i)+twopi
    end if
! calculate the piezoelectric tensor element from difference in polarisations
    t1=wkptnr*occmax*dble(nkspolar*ngridk(i))/(twopi*deltast)
    pelt(i,istrain)=t1*(pvl2(i)-pvl1(i))
  end do
end do
if (mp_mpi) then
  open(50,file='PIEZOELT.OUT',form='FORMATTED')
  write(50,*)
  write(50,'("Lattice vector matrix, A, changed by")')
  write(50,*)
  write(50,'("     A → A + eₖ dt,")')
  write(50,*)
  write(50,'("where dt is an infinitesimal scalar and eₖ is a strain tensor")')
  write(50,*)
  write(50,'("The piezoelectric tensor is the derivative of the polarisation &
   &vector dPᵢ/dt, i=1...3")')
  do k=1,nstrain
    write(50,*)
    write(50,'("Strain tensor k : ",I1)') k
    do j=1,3
      write(50,'(3G18.10)') (strain(i,j,k),i=1,3)
    end do
    write(50,'("Piezoelectric tensor components dPᵢ/dt, i=1...3 :")')
    write(50,'(" lattice coordinates : ",3G18.10)') pelt(:,k)
    call r3mv(avec,pelt(:,k),vc)
    write(50,'(" Cartesian coordinates : ",3G18.10)') vc(:)
    write(50,'("  length : ",G18.10)') norm2(vc(:))
  end do
  close(50)
  write(*,*)
  write(*,'("Info(piezoelt):")')
  write(*,'(" Piezoelectric tensor written to PIEZOELT.OUT")')
end if
! write test file if required
call writetest(380,'Piezoelectric tensor',nv=3*nstrain,tol=1.d-5,rva=pelt)
deallocate(pelt)
! restore original parameters
tshift=tshift0
avec(:,:)=avec0(:,:)
end subroutine

