
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: sort
! !INTERFACE:
subroutine sort(n,x)
! !INPUT/OUTPUT PARAMETERS:
!   n : number of elements in array (in,integer)
!   x : real array (inout,real(n))
! !DESCRIPTION:
!   Sorts elements of a real array into ascending order. See {\tt sortidx}.
!
! !REVISION HISTORY:
!   Created May 2024 (JKD)
!EOP
!BOC
implicit none
integer, intent(in) :: n
real(8), intent(inout) :: x(n)
! automatic arrays
integer idx(n)
real(8) xt(n)
xt(1:n)=x(1:n)
! find the permutation index
call sortidx(n,xt,idx)
x(1:n)=xt(idx(1:n))
end subroutine
!EOC

