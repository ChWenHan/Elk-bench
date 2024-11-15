
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: occupy
! !INTERFACE:
subroutine occupy
! !USES:
use modmain
use modtest
! !DESCRIPTION:
!   Finds the Fermi energy and sets the occupation numbers for the
!   second-variational states using the routine {\tt fermi}.
!
! !REVISION HISTORY:
!   Created February 2004 (JKD)
!   Added gap estimation, November 2009 (F. Cricchio)
!   Added adaptive smearing width, April 2010 (T. Bjorkman)
!EOP
!BOC
implicit none
! local variables
integer, parameter :: maxit=1000
integer ik,ist,it
real(8) e0,e1,e,de
real(8) ei0,ei1,ed0,ed1
real(8) chg,w,x,t0
! external functions
real(8), external :: sdelta,stheta
! determine the smearing width automatically if required
if ((autoswidth).and.(iscl > 1)) call findswidth
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
if (e0 < e0min) then
  write(*,*)
  write(*,'("Warning(occupy): minimum eigenvalue less than minimum &
   &linearisation energy : ",2G18.10)') e0,e0min
  write(*,'(" for s.c. loop ",I5)') iscl
end if
t0=1.d0/swidth
! determine the Fermi energy using the bisection method
do it=1,maxit
  efermi=0.5d0*(e0+e1)
  chg=0.d0
  do ik=1,nkpt
    w=wkpt(ik)
    do ist=1,nstsv
      e=evalsv(ist,ik)
      if (e < e0min) then
        occsv(ist,ik)=0.d0
      else
        x=(efermi-e)*t0
        occsv(ist,ik)=occmax*stheta(stype,x)
        chg=chg+w*occsv(ist,ik)
      end if
    end do
  end do
  if (chg < chgval) then
    e0=efermi
  else
    e1=efermi
  end if
  if ((e1-e0) < 1.d-12) goto 10
end do
write(*,*)
write(*,'("Warning(occupy): could not find Fermi energy")')
10 continue
! find the density of states at the Fermi surface in units of
! states/Hartree/unit cell
fermidos=0.d0
do ik=1,nkpt
  w=wkpt(ik)
  do ist=1,nstsv
    x=(evalsv(ist,ik)-efermi)*t0
    fermidos=fermidos+w*sdelta(stype,x)
  end do
  if (abs(occsv(nstsv,ik)) > epsocc) then
    write(*,*)
    write(*,'("Warning(occupy): not enough empty states for k-point ",I6)') ik
    write(*,'(" and s.c. loop ",I5)') iscl
  end if
end do
fermidos=fermidos*occmax*t0
! write Fermi density of states to test file
call writetest(500,'DOS at Fermi energy',tol=5.d-3,rv=fermidos)
! estimate the indirect and direct band gaps (FC)
ei0=-1.d8; ei1=1.d8
de=1.d8
ikgap(1:3)=1
do ik=1,nkpt
  ed0=-1.d8; ed1=1.d8
  do ist=1,nstsv
    e=evalsv(ist,ik)
    if (e <= efermi) then
      if (e > ed0) ed0=e
      if (e > ei0) then
! transfer is a workaround for a bug in Intel Fortran versions 17 and 18
        ikgap(1)=transfer(ik,ik)
        ei0=e
      end if
    else
      if (e < ed1) ed1=e
      if (e < ei1) then
        ikgap(2)=ik
        ei1=e
      end if
    end if
  end do
  e=ed1-ed0
  if (e < de) then
    ikgap(3)=ik
    de=e
  end if
end do
bandgap(1)=ei1-ei0
bandgap(2)=de
! write band gap to test file
call writetest(510,'estimated indirect band gap',tol=2.d-2,rv=bandgap(1))
end subroutine
!EOC

