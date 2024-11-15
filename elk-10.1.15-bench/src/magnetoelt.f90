
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine magnetoelt
use modmain
use modmpi
use modtest
implicit none
! local variables
integer i,j
real(8) pvl1(3),pvl2(3)
real(8) felt(3,3),vc(3),t1
! initialise universal variables
call init0
call init1
! store original parameters
bfieldc00(:)=bfieldc0(:)
ngridk0(:)=ngridk(:)
maxscl0=maxscl
spinpol0=spinpol
reducebf0=reducebf
tshift0=tshift
! no shifting of the atomic basis
tshift=.false.
! enable spin-polarisation
spinpol=.true.
! no magnetic field reduction
reducebf=1.d0
! initial ground-state run should start from atomic densities
trdstate=.false.
do j=1,3
  if (mp_mpi) then
    write(*,'("Info(magnetoelt): working on magnetic field component ",I1)') j
  end if
! apply negative change to magnetic field
  bfieldc0(j)=bfieldc00(j)-0.5d0*deltabf
! run the ground-state calculation
  call gndstate
! subsequent calculations will read in the previous potential
  trdstate=.true.
! compute the first polarisation in lattice coordinates
  call polar(pvl1)
! apply positive change to magnetic field
  bfieldc0(j)=bfieldc00(j)+0.5d0*deltabf
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
! calculate the magnetoelectric tensor element from difference in polarisations
    t1=wkptnr*occmax*dble(nkspolar*ngridk(i))/(twopi*deltabf)
    felt(i,j)=t1*(pvl2(i)-pvl1(i))
  end do
end do
if (mp_mpi) then
  open(50,file='MAGNETOELT.OUT',form='FORMATTED')
  write(50,*)
  write(50,'("The magnetoelectric tensor is the change in the polarisation")')
  write(50,'("with respect to the external magnetic field dPᵢ/dBⱼ, for")')
  write(50,'("components i,j=1...3")')
  do j=1,3
    write(50,*)
    write(50,'("Magnetic field Cartesian component j : ",I1)') j
    write(50,'("Magnetoelectric tensor components dPᵢ/dBⱼ, i=1...3")')
    write(50,'(" lattice coordinates : ",3G18.10)') felt(:,j)
    call r3mv(avec,felt(:,j),vc)
    write(50,'(" Cartesian coordinates : ",3G18.10)') vc(:)
    write(50,'("  length : ",G18.10)') norm2(vc(:))
  end do
  close(50)
  write(*,*)
  write(*,'("Info(magnetoelt):")')
  write(*,'(" Magnetoelectric tensor writtent to MAGNETOELT.OUT")')
end if
! write test file if required
call writetest(390,'Magnetoelectric tensor',nv=9,tol=1.d-5,rva=felt)
! restore original parameters
tshift=tshift0
spinpol=spinpol0
bfieldc0(:)=bfieldc00(:)
reducebf=reducebf0
end subroutine

