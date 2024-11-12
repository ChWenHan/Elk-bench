
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dhmlaa(is,ias,ngp,ngpq,apwalm,apwalmq,dapwalm,dapwalmq,ld,dh)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: is,ias
integer, intent(in) :: ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: dh(*)
! local variables
integer io,jo,i
integer l0,l1,l2,l3
integer lm1,lm3,lma,lmb
real(8) t0,t1
complex(8) z1
! automatic arrays
complex(8) y1(ngp),a1(lmoapw(is),ngpq),b1(lmoapw(is),ngp)
complex(8) y2(ngp),a2(lmoapw(is),ngpq),b2(lmoapw(is),ngp)
t0=0.5d0*rmt(is)**2
i=0
do l1=0,lmaxapw
  do lm1=l1**2+1,(l1+1)**2
    do io=1,apword(l1,is)
      i=i+1
      t1=t0*apwfr(nrmt(is),1,io,l1,ias)
      y1(1:ngp)=0.d0
      do l3=0,lmaxapw
        if (mod(l1+l3,2) == 0) then; l0=0; else; l0=1; end if
        do lm3=l3**2+1,(l3+1)**2
          do jo=1,apword(l3,is)
            z1=0.d0
            do l2=l0,lmaxo,2
              lma=l2**2+1; lmb=lma+2*l2
              z1=z1+sum(gntyyy(lma:lmb,lm3,lm1)*dhaa(lma:lmb,jo,l3,io,l1,ias))
            end do
            if (abs(dble(z1))+abs(aimag(z1)) > 1.d-14) then
              call zaxpy(ngp,z1,apwalm(:,jo,lm3),1,y1,1)
            end if
          end do
        end do
      end do
      if (ias == iasph) then
        y2(1:ngp)=0.d0
        do l3=0,lmaxapw
          if (mod(l1+l3,2) == 0) then; l0=0; else; l0=1; end if
          do lm3=l3**2+1,(l3+1)**2
            do jo=1,apword(l3,is)
              z1=0.d0
              do l2=l0,lmaxo,2
                lma=l2**2+1; lmb=lma+2*l2
                z1=z1+sum(gntyry(lma:lmb,lm3,lm1)*haa(lma:lmb,jo,l3,io,l1,ias))
              end do
! kinetic surface contribution
              if (lm1 == lm3) z1=z1+t1*apwdfr(jo,l1,ias)
              if (abs(dble(z1))+abs(aimag(z1)) > 1.d-14) then
                call zaxpy(ngp,z1,dapwalm(:,jo,lm3),1,y1,1)
                call zaxpy(ngp,z1,apwalm(:,jo,lm3),1,y2,1)
              end if
            end do
          end do
        end do
        a2(i,1:ngpq)=dapwalmq(1:ngpq,io,lm1)
        b2(i,1:ngp)=y2(1:ngp)
      end if
      a1(i,1:ngpq)=apwalmq(1:ngpq,io,lm1)
      b1(i,1:ngp)=y1(1:ngp)
    end do
  end do
end do
call zmctm(lmoapw(is),ngpq,ngp,a1,b1,ld,dh)
if (ias == iasph) then
  call zmctm(lmoapw(is),ngpq,ngp,a2,b2,ld,dh)
end if
end subroutine

