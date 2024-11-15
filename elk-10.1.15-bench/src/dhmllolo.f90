
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine dhmllolo(is,ias,ngp,ngpq,ld,dh)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: is,ias
integer, intent(in) :: ngp,ngpq
integer, intent(in) :: ld
complex(8), intent(out) :: dh(ld,*)
! local variables
integer ilo,jlo,i,j
integer l0,l1,l2,l3
integer lm1,lm3,lma,lmb
complex(8) z1
do jlo=1,nlorb(is)
  l3=lorbl(jlo,is)
  do lm3=l3**2+1,(l3+1)**2
    j=ngp+idxlo(lm3,jlo,ias)
    do ilo=1,nlorb(is)
      l1=lorbl(ilo,is)
      if (mod(l1+l3,2) == 0) then; l0=0; else; l0=1; end if
      do lm1=l1**2+1,(l1+1)**2
        i=ngpq+idxlo(lm1,ilo,ias)
        z1=0.d0
        do l2=l0,lmaxo,2
          lma=l2**2+1; lmb=(l2+1)**2
          z1=z1+sum(gntyyy(lma:lmb,lm3,lm1)*dhlolo(lma:lmb,jlo,ilo,ias))
        end do
        dh(i,j)=z1
      end do
    end do
  end do
end do
end subroutine

