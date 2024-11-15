
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zlrzncnv(n,s,w,zf)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: s,w(n)
complex(8), intent(inout) :: zf(n)
! local variables
integer i,j
real(8), parameter :: pi=3.1415926535897932385d0
real(8) s2,wi,dw
complex(8) zsm
! automatic arrays
complex(8) zg(n)
s2=s**2
do i=1,n
  wi=w(i)
  zsm=0.d0
  do j=1,n-1
    dw=w(j+1)-w(j)
    zsm=zsm+dw*zf(j)/((w(j)-wi)**2+s2)
  end do
  zg(i)=zsm*s/pi
end do
zf(1:n)=zg(1:n)
end subroutine

