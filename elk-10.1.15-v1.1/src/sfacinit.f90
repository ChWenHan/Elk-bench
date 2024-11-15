
! Copyright (C) 2010 A. I. Baranov and F. Wagner.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine sfacinit
use modmain
use modpw
implicit none
! local variables
logical trhonorm0
integer ik,ist,is,ias
! allocatable arrays
real(8), allocatable :: occcr0(:,:)
! initialise universal variables
call init0
call init1
! read density and potentials from file
call readstate
! use existing density if wsfac is default
if ((wsfac(1) <= -1.d6).and.(wsfac(2) >= 1.d6)) goto 10
! make a copy of the core state occupation numbers
allocate(occcr0(nstspmax,natmtot))
occcr0(:,:)=occcr(:,:)
! zero the core state occupation numbers for eigenvalues not in energy window
do ias=1,natmtot
  is=idxis(ias)
  do ist=1,nstsp(is)
    if (spcore(ist,is)) then
      if ((evalcr(ist,ias) < wsfac(1)).or.(evalcr(ist,ias) > wsfac(2))) then
        occcr(ist,ias)=0.d0
      end if
    end if
  end do
end do
! generate the core wavefunctions and densities
call gencore
! restore the core state occupation numbers
occcr(:,:)=occcr0(:,:)
deallocate(occcr0)
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
do ik=1,nkpt
! get the eigenvalues and occupation numbers from file
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
! zero occupation numbers for eigenvalues not in energy window
  do ist=1,nstsv
    if ((evalsv(ist,ik) < wsfac(1)).or.(evalsv(ist,ik) > wsfac(2))) then
      occsv(ist,ik)=0.d0
    end if
  end do
end do
! computed density should not be normalised
trhonorm0=trhonorm
trhonorm=.false.
! generate the density and magnetisation
call rhomag
trhonorm=trhonorm0
10 continue
! generate the H-vectors
call genhvec
end subroutine

