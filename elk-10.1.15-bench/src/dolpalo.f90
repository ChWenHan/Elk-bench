
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine dolpalo(is,ias,ngp,ngpq,dapwalm,dapwalmq,ld,od)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: is,ias,ngp,ngpq
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: od(ld,*)
! local variables
integer ilo,io,l,lm,i,j,k
real(8) t1
if (ias /= iasph) return
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do lm=l**2+1,(l+1)**2
    k=idxlo(lm,ilo,ias)
    i=ngpq+k
    j=ngp+k
    do io=1,apword(l,is)
      t1=oalo(io,ilo,ias)
      od(1:ngpq,j)=od(1:ngpq,j)+t1*conjg(dapwalmq(1:ngpq,io,lm))
      od(i,1:ngp)=od(i,1:ngp)+t1*dapwalm(1:ngp,io,lm)
    end do
  end do
end do
end subroutine

