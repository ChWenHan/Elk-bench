
! Copyright (C) 2022 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rzmctmu(l,n,a,b,ld,c)
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
call holdthd(n,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i) &
!$OMP NUM_THREADS(nthd)
do j=1,n
  do i=1,j
    c(i,j)=c(i,j)+dot_product(a(:,i),b(:,j))
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

