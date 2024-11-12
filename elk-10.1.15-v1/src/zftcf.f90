
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zftcf(ngp,jlgpr,ylmgp,ld,sfacgp,cfmt,cfir,zfgp)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: jlgpr(njcmax,nspecies,ngp)
complex(8), intent(in) :: ylmgp(lmmaxo,ngp)
integer, intent(in) :: ld
complex(8), intent(in) :: sfacgp(ld,natmtot)
complex(4), intent(in) :: cfmt(npcmtmax,natmtot),cfir(ngtc)
complex(8), intent(out) :: zfgp(ngp)
! local variables
integer is,ia,ias,ig
integer nrc,nrci,irco,irc
integer l,lm,n,i,j
real(8) t0,y0
complex(8) zsm,z1
! automatic arrays
complex(4) ylm(lmmaxo),cfft(ngtc)
!-----------------------------------!
!     interstitial contribution     !
!-----------------------------------!
! multiply by coarse characteristic function
cfft(1:ngtc)=cfir(1:ngtc)*cfrc(1:ngtc)
! Fourier transform to coarse G-grid
call cfftifc(3,ngdgc,-1,cfft)
zfgp(1:ngp)=cfft(igfc(1:ngp))
!---------------------------------!
!     muffin-tin contribution     !
!---------------------------------!
t0=1.d0/omega
y0=t0*fourpi*y00
do ig=1,ngp
  ylm(1:lmmaxo)=t0*ylmgp(1:lmmaxo,ig)
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    irco=nrci+1
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      zsm=0.d0
      i=1
      j=1
! inner part of muffin-tin, note that lmaxi >= 1
      if (lmaxi == 1) then
        do irc=1,nrci
          zsm=zsm+wrcmt(irc,is)*(jlgpr(j,is,ig)*cfmt(i,ias)*y0+jlgpr(j+1,is,ig)&
           *(cfmt(i+1,ias)*ylm(2)+cfmt(i+2,ias)*ylm(3)+cfmt(i+3,ias)*ylm(4)))
          i=i+4
          j=j+2
        end do
      else
        do irc=1,nrci
          z1=jlgpr(j,is,ig)*cfmt(i,ias)*y0+jlgpr(j+1,is,ig) &
           *(cfmt(i+1,ias)*ylm(2)+cfmt(i+2,ias)*ylm(3)+cfmt(i+3,ias)*ylm(4))
          i=i+4
          j=j+2
          do l=2,lmaxi
            n=2*l
            lm=l**2+1
            z1=z1+jlgpr(j,is,ig)*sum(cfmt(i:i+n,ias)*ylm(lm:lm+n))
            i=i+n+1
            j=j+1
          end do
          zsm=zsm+wrcmt(irc,is)*z1
        end do
      end if
! outer part of muffin-tin, note that lmaxo >= 3
      do irc=irco,nrc
        z1=jlgpr(j,is,ig)*cfmt(i,ias)*y0+jlgpr(j+1,is,ig) &
         *(cfmt(i+1,ias)*ylm(2)+cfmt(i+2,ias)*ylm(3)+cfmt(i+3,ias)*ylm(4)) &
         +jlgpr(j+2,is,ig) &
         *(cfmt(i+4,ias)*ylm(5)+cfmt(i+5,ias)*ylm(6)+cfmt(i+6,ias)*ylm(7) &
          +cfmt(i+7,ias)*ylm(8)+cfmt(i+8,ias)*ylm(9)) &
         +jlgpr(j+3,is,ig) &
         *(cfmt(i+9,ias)*ylm(10)+cfmt(i+10,ias)*ylm(11)+cfmt(i+11,ias)*ylm(12) &
          +cfmt(i+12,ias)*ylm(13)+cfmt(i+13,ias)*ylm(14)+cfmt(i+14,ias)*ylm(15)&
          +cfmt(i+15,ias)*ylm(16))
        i=i+16
        j=j+4
        do l=4,lmaxo
          n=2*l
          lm=l**2+1
          z1=z1+jlgpr(j,is,ig)*sum(cfmt(i:i+n,ias)*ylm(lm:lm+n))
          i=i+n+1
          j=j+1
        end do
        zsm=zsm+wrcmt(irc,is)*z1
      end do
      zfgp(ig)=zfgp(ig)+conjg(sfacgp(ig,ias))*zsm
    end do
  end do
end do
end subroutine

