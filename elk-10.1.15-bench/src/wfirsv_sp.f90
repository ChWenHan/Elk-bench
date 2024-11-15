
! Copyright (C) 2022 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wfirsv_sp(tgp,nst,idx,ngdg,igf,ngp,igpig,evecfv,evecsv,ld,wfir)
use modmain
use modomp
implicit none
! arguments
logical, intent(in) :: tgp
integer, intent(in) :: nst,idx(*),ngdg(3),igf(*)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
integer, intent(in) :: ld
complex(4), intent(out) :: wfir(ld,nspinor,nst)
! local variables
integer ist,ispn,jspn
integer n,igp,i,j,k,nthd
real(8) t0
complex(8) z1
! automatic arrays
complex(4) wfgp(ngkmax)
t0=1.d0/sqrt(omega)
call holdthd(nst,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(wfgp,k,ispn,jspn) &
!$OMP PRIVATE(n,ist,i,z1,igp) &
!$OMP NUM_THREADS(nthd)
do j=1,nst
! index to state in evecsv
  if (idx(1) == 0) then; k=j; else; k=idx(j); end if
  if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
    do ispn=1,nspinor
      jspn=jspnfv(ispn)
      n=ngp(jspn)
      if (tgp) then
        wfir(1:n,ispn,j)=0.e0
      else
        wfgp(1:n)=0.e0
      end if
      do ist=1,nstfv
        i=(ispn-1)*nstfv+ist
        z1=evecsv(i,k)
        if (abs(dble(z1))+abs(aimag(z1)) < epsocc) cycle
        if (tgp) then
! wavefunction in G+p-space
          wfir(1:n,ispn,j)=wfir(1:n,ispn,j)+z1*evecfv(1:n,ist,jspn)
        else
! wavefunction in real-space
          z1=t0*z1
          wfgp(1:n)=wfgp(1:n)+z1*evecfv(1:n,ist,jspn)
        end if
      end do
! Fourier transform wavefunction to real-space if required
      if (.not.tgp) then
        wfir(1:ld,ispn,j)=0.e0
        do igp=1,n
          wfir(igf(igpig(igp,jspn)),ispn,j)=wfgp(igp)
        end do
        call cfftifc(3,ngdg,1,wfir(:,ispn,j))
      end if
    end do
  else
! spin-unpolarised wavefunction
    n=ngp(1)
    if (tgp) then
      wfir(1:n,1,j)=evecfv(1:n,k,1)
    else
      wfir(1:ld,1,j)=0.e0
      do igp=1,n
        wfir(igf(igpig(igp,1)),1,j)=t0*evecfv(igp,k,1)
      end do
      call cfftifc(3,ngdg,1,wfir(:,1,j))
    end if
  end if
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

