
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine gencrm(n,wf11,wf12,wf21,wf22,crho,ld,cmag)
use modmain
implicit none
! arguments
integer, intent(in) :: n
complex(4), intent(in) :: wf11(n),wf12(n),wf21(n),wf22(n)
complex(4), intent(out) :: crho(n)
integer, intent(in) :: ld
complex(4), intent(out) :: cmag(ld,ndmag)
! local variables
integer i
complex(4) c11,c12,c21,c22,c1,c2
if (ncmag) then
! non-collinear case
!$OMP SIMD PRIVATE(c11,c12,c21,c22,c1,c2) SIMDLEN(8)
  do i=1,n
    c11=wf11(i); c12=wf12(i)
    c21=wf21(i); c22=wf22(i)
! up-dn spin density
    c1=conjg(c11)*c22
! dn-up spin density
    c2=conjg(c12)*c21
! x-component: up-dn + dn-up
    cmag(i,1)=c1+c2
! y-component: i*(dn-up - up-dn)
    c1=c2-c1
    cmag(i,2)=cmplx(-aimag(c1),real(c1),4)
    c1=conjg(c11)*c21
    c2=conjg(c12)*c22
! z-component: up-up - dn-dn
    cmag(i,3)=c1-c2
! density: up-up + dn-dn
    crho(i)=c1+c2
  end do
else
! collinear case
!$OMP SIMD PRIVATE(c1,c2) SIMDLEN(8)
  do i=1,n
    c1=conjg(wf11(i))*wf21(i)
    c2=conjg(wf12(i))*wf22(i)
    cmag(i,1)=c1-c2
    crho(i)=c1+c2
  end do
end if
end subroutine

