
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine olpalo(is,ias,ngp,apwalm,ld,o)
use modmain
implicit none
! arguments
integer, intent(in) :: is,ias,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: o(ld,*)
! local variables
integer ilo,io,l,lm,j
real(8) t1
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do lm=l**2+1,(l+1)**2
    j=ngp+idxlo(lm,ilo,ias)
    do io=1,apword(l,is)
      t1=oalo(io,ilo,ias)
      o(1:ngp,j)=o(1:ngp,j)+t1*conjg(apwalm(1:ngp,io,lm))
    end do
  end do
end do
end subroutine

