
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlaaq(is,ias,ngp,ngpq,apwalm,apwalmq,ld,hq)
use modmain
implicit none
integer, intent(in) :: is,ias,ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: hq(*)
! local variables
integer io,jo,i
integer l0,l1,l2,l3
integer lm1,lm3,lma,lmb
real(8) t0,t1
complex(8) z1
! automatic arrays
complex(8) y(ngp),a(lmoapw(is),ngpq),b(lmoapw(is),ngp)
t0=0.5d0*rmt(is)**2
i=0
do l1=0,lmaxapw
  do lm1=l1**2+1,(l1+1)**2
    do io=1,apword(l1,is)
      i=i+1
      t1=t0*apwfr(nrmt(is),1,io,l1,ias)
      y(:)=0.d0
      do l3=0,lmaxapw
        if (mod(l1+l3,2) == 0) then; l0=0; else; l0=1; end if
        do lm3=l3**2+1,(l3+1)**2
          do jo=1,apword(l3,is)
            z1=0.d0
! kinetic and potential contribution
            do l2=l0,lmaxo,2
              lma=l2**2+1; lmb=lma+2*l2
              z1=z1+sum(gntyry(lma:lmb,lm3,lm1)*haa(lma:lmb,jo,l3,io,l1,ias))
            end do
! kinetic surface contribution
            if (lm1 == lm3) z1=z1+t1*apwdfr(jo,l1,ias)
            if (abs(dble(z1))+abs(aimag(z1)) > 1.d-14) then
              call zaxpy(ngp,z1,apwalm(:,jo,lm3),1,y,1)
            end if
          end do
        end do
      end do
      a(i,1:ngpq)=apwalmq(1:ngpq,io,lm1)
      b(i,1:ngp)=y(1:ngp)
    end do
  end do
end do
call zmctm(lmoapw(is),ngpq,ngp,a,b,ld,hq)
end subroutine

