
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine cftwfir(ngp,igpig,wfir)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(4), intent(inout) :: wfir(ngtc,nspinor,nstsv)
! local variables
integer ist,ispn,jspn
integer n,igp,nthd
real(8) t0
! automatic arrays
complex(4) c(ngkmax)
t0=1.d0/sqrt(omega)
call holdthd(nstsv,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(c,ispn,jspn,n,igp) &
!$OMP NUM_THREADS(nthd)
do ist=1,nstsv
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    n=ngp(jspn)
    c(1:n)=wfir(1:n,ispn,ist)
    wfir(1:ngtc,ispn,ist)=0.e0
    do igp=1,n
      wfir(igfc(igpig(igp,jspn)),ispn,ist)=t0*c(igp)
    end do
    call cfftifc(3,ngdgc,1,wfir(:,ispn,ist))
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

