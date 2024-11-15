
! Copyright (C) 2005 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: fsmooth
! !INTERFACE:
pure subroutine fsmooth(m,n,f)
! !INPUT/OUTPUT PARAMETERS:
!   m  : number of 3-point running averages to perform (in,integer)
!   n  : number of point (in,integer)
!   f  : function array (inout,real(n))
! !DESCRIPTION:
!   Removes numerical noise from a function by performing $m$ successive
!   3-point running averages on the data. The endpoints are kept fixed.
!
! !REVISION HISTORY:
!   Created December 2005 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m,n
real(8), intent(inout) :: f(n)
! local variables
integer i,j
real(8) f1,f2,f3
do i=1,m
  f1=f(1); f2=f(2)
  do j=2,n-1
    f3=f(j+1)
    f(j)=0.3333333333333333333d0*(f1+f2+f3)
    f1=f2; f2=f3
  end do
end do
end subroutine
!EOC

