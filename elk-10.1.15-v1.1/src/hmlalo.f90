
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine hmlalo(is,ias,ngp,apwalm,ld,h)
use modmain
implicit none
! arguments
integer, intent(in) :: is,ias,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: h(ld,*)
! local variables
integer io,ilo,j
integer l0,l1,l2,l3
integer lm1,lm3,lma,lmb
complex(8) z1
do ilo=1,nlorb(is)
  l1=lorbl(ilo,is)
  do lm1=l1**2+1,(l1+1)**2
    j=ngp+idxlo(lm1,ilo,ias)
    do l3=0,lmaxapw
      if (mod(l1+l3,2) == 0) then; l0=0; else; l0=1; end if
      do lm3=l3**2+1,(l3+1)**2
        do io=1,apword(l3,is)
          z1=0.d0
          do l2=l0,lmaxo,2
            lma=l2**2+1; lmb=lma+2*l2
            z1=z1+sum(gntyry(lma:lmb,lm3,lm1)*hloa(lma:lmb,io,l3,ilo,ias))
          end do
! note that what is actually computed is the Hermitian conjugate of <lo|H|APW>
          if (abs(dble(z1))+abs(aimag(z1)) > 1.d-14) then
            h(1:ngp,j)=h(1:ngp,j)+conjg(z1*apwalm(1:ngp,io,lm3))
          end if
        end do
      end do
    end do
  end do
end do
end subroutine

