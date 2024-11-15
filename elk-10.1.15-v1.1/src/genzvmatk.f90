
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzvmatk(zvmt,zvir,ngp,igpig,wfmt,wfir,wfgp,vmat)
use modmain
use modomp
implicit none
! arguments
! the potential is multiplied by the radial integration weights in the
! muffin-tin and by the characteristic function in the interstitial region
complex(8), intent(in) :: zvmt(npcmtmax,natmtot),zvir(ngtc)
integer, intent(in) :: ngp,igpig(ngp)
complex(4), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
! note that wfir does not have a 1/sqrt(omega) prefactor
complex(4), intent(in) :: wfir(ngtc,nspinor,nstsv)
complex(4), intent(in) :: wfgp(ngp,nspinor,nstsv)
complex(8), intent(out) :: vmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,nthd
integer is,ias,npc,igp
! automatic arrays
complex(4) wfmt1(npcmtmax),wfir1(ngtc),c(ngp)
! external functions
complex(4), external :: cdotc
call holdthd(nstsv,nthd)
! zero the matrix elements
vmat(1:nstsv,1:nstsv)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,ispn,ias) &
!$OMP PRIVATE(is,npc,ist) &
!$OMP NUM_THREADS(nthd)
do jst=1,nstsv
  do ispn=1,nspinor
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
! apply complex potential to wavefunction
      wfmt1(1:npc)=zvmt(1:npc,ias)*wfmt(1:npc,ias,ispn,jst)
! compute the inner products
      do ist=1,nstsv
        vmat(ist,jst)=vmat(ist,jst)+cdotc(npc,wfmt(:,ias,ispn,ist),1,wfmt1,1)
      end do
    end do
  end do
end do
!$OMP END PARALLEL DO
!---------------------------!
!     interstitial part     !
!---------------------------!
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,c,ispn,igp,ist) &
!$OMP NUM_THREADS(nthd)
do jst=1,nstsv
  do ispn=1,nspinor
! apply potential to wavefunction
    wfir1(1:ngtc)=zvir(1:ngtc)*wfir(1:ngtc,ispn,jst)
! Fourier transform to G+p-space
    call cfftifc(3,ngdgc,-1,wfir1)
    do igp=1,ngp
      c(igp)=wfir1(igfc(igpig(igp)))
    end do
    do ist=1,nstsv
! compute inner product
      vmat(ist,jst)=vmat(ist,jst)+cdotc(ngp,wfgp(:,ispn,ist),1,c,1)
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

