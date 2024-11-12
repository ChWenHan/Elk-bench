
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rfinterp
! !INTERFACE:
subroutine rfinterp(ni,xi,wci,fi,no,xo,fo)
! !INPUT/OUTPUT PARAMETERS:
!   ni  : number of input points (in,integer)
!   xi  : input abscissa array (in,real(ni))
!   wci : input spline coefficient weights (in,real(12,ni))
!   fi  : input data array (in,real(ni))
!   no  : number of output points (in,integer)
!   xo  : output abscissa array (in,real(no))
!   fo  : output interpolated function (out,real(no))
! !DESCRIPTION:
!   Given a function defined on a set of input points, this routine uses a
!   clamped cubic spline to interpolate the function on a different set of
!   points. See routine {\tt spline}.
!
! !REVISION HISTORY:
!   Created January 2005 (JKD)
!   Arguments changed, April 2016 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ni
real(8), intent(in) :: xi(ni),wci(12,ni),fi(ni)
integer, intent(in) :: no
real(8), intent(in) :: xo(no)
real(8), intent(out) :: fo(no)
! local variables
integer i,j,k,l
real(8) x,dx
! automatic arrays
real(8) cf(3,ni)
if (ni == 1) then
  fo(1:no)=fi(1)
  return
end if
! compute the spline coefficients
call splinew(ni,wci,fi,cf)
! evaluate spline at output points
i=1
do l=1,no
  x=xo(l)
  if (i >= ni) i=1
  if (x >= xi(i)) then
    if (x > xi(i+1)) then
! binary search
      i=1
      j=ni
      do while (j > i+1)
        k=(i+j)/2
        if (x < xi(k)) then
          j=k
        else
          i=k
        end if
      end do
    end if
  end if
  dx=x-xi(i)
  fo(l)=fi(i)+dx*(cf(1,i)+dx*(cf(2,i)+dx*cf(3,i)))
end do
end subroutine
!EOC

