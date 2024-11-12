
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zmctmu(l,n,a,b,ld,c)
use modomp
implicit none
! arguments
integer, intent(in) :: l,n
! pass in complex arrays as real
real(8), intent(in) :: a(2*l,n),b(2*l,n)
integer, intent(in) :: ld
complex(8), intent(inout) :: c(ld,*)
! local variables
integer i,j,nthd
! external functions
complex(8), external :: zdotc
call holdthd(n,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i) &
!$OMP NUM_THREADS(nthd)
do j=1,n
  do i=1,j-1
    c(i,j)=c(i,j)+zdotc(l,a(1,i),1,b(1,j),1)
  end do
  c(j,j)=c(j,j)+dot_product(a(:,j),b(:,j))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

