
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: sortidx
! !INTERFACE:
pure subroutine sortidx(n,x,idx)
! !INPUT/OUTPUT PARAMETERS:
!   n   : number of elements in array (in,integer)
!   x   : real array (in,real(n))
!   idx : permutation index (out,integer(n))
! !DESCRIPTION:
!   Finds the permutation index {\tt idx} which sorts the real array {\tt x}
!   into ascending order. No sorting of the array {\tt x} itself is performed.
!   Uses the heapsort algorithm.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!   Included tolerance eps, April 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n)
integer, intent(out) :: idx(n)
! local variables
integer i,j,k,l,m
! tolerance for deciding if one number is smaller than another
real(8), parameter :: eps=1.d-14
if (n <= 0) return
do i=1,n
  idx(i)=i
end do
if (n == 1) return
l=n/2+1
k=n
do
  if (l > 1) then
    l=l-1
    m=idx(l)
  else
    m=idx(k)
    idx(k)=idx(1)
    k=k-1
    if (k == 1) then
      idx(1)=m
      return
    end if
  end if
  i=l
  j=l+l
  do while (j <= k)
    if (j < k) then
      if (x(idx(j)) < x(idx(j+1))+eps) j=j+1
    end if
    if (x(m) < x(idx(j))+eps) then
      idx(i)=idx(j)
      i=j
      j=j+j
    else
      j=k+1
    end if
  end do
  idx(i)=m
end do
end subroutine
!EOC

