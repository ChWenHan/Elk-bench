
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlaloq(is,ias,ngp,ngpq,apwalm,apwalmq,ld,hq)
use modmain
implicit none
! arguments
integer, intent(in) :: is,ias
integer, intent(in) :: ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: hq(ld,*)
! local variables
integer io,ilo
integer l0,l1,l2,l3
integer lm1,lm3,lma,lmb
integer i0,j0,i,j
complex(8) z1
do ilo=1,nlorb(is)
  l1=lorbl(ilo,is)
  do lm1=l1**2+1,(l1+1)**2
    i=idxlo(lm1,ilo,ias)
    i0=ngpq+i
    j0=ngp+i
    do l3=0,lmaxapw
      if (mod(l1+l3,2) == 0) then; l0=0; else; l0=1; end if
      do lm3=l3**2+1,(l3+1)**2
        do io=1,apword(l3,is)
          z1=0.d0
          do l2=l0,lmaxo,2
            lma=l2**2+1; lmb=lma+2*l2
            z1=z1+sum(gntyry(lma:lmb,lm3,lm1)*hloa(lma:lmb,io,l3,ilo,ias))
          end do
          if (abs(dble(z1))+abs(aimag(z1)) > 1.d-14) then
            do i=1,ngpq
              hq(i,j0)=hq(i,j0)+conjg(z1*apwalmq(i,io,lm3))
            end do
            do j=1,ngp
              hq(i0,j)=hq(i0,j)+z1*apwalm(j,io,lm3)
            end do
          end if
        end do
      end do
    end do
  end do
end do
end subroutine

