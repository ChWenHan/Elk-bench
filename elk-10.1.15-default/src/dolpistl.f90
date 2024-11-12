
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine dolpistl(ngp,ngpq,igpig,igpqig,ld,od)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ngp,ngpq
integer, intent(in) :: igpig(ngkmax),igpqig(ngkmax)
integer, intent(in) :: ld
complex(8), intent(out) :: od(ld,*)
! local variables
integer ig,j1,j2,j3,i,j
do j=1,ngp
  ig=igpig(j)
  j1=ivg(1,ig); j2=ivg(2,ig); j3=ivg(3,ig)
  do i=1,ngpq
    ig=igpqig(i)
    ig=ivgig(ivg(1,ig)-j1,ivg(2,ig)-j2,ivg(3,ig)-j3)
    od(i,j)=dcfunig(ig)
  end do
end do
end subroutine

