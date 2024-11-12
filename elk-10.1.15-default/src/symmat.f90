
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symmat(al)
use modmain
implicit none
! arguments
real(8), intent(inout) :: al(3,3)
! local variables
integer isym,lspl
real(8) as(3,3),s(3,3)
real(8) b(3,3),c(3,3),t1
! make symmetric average
as(1:3,1:3)=0.d0
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  s(1:3,1:3)=dble(symlat(1:3,1:3,lspl))
  call r3mtm(s,al,b)
  call r3mm(b,s,c)
  as(1:3,1:3)=as(1:3,1:3)+c(1:3,1:3)
end do
! normalise
t1=1.d0/dble(nsymcrys)
al(1:3,1:3)=t1*as(1:3,1:3)
end subroutine

