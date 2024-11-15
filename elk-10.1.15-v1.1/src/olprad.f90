
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: olprad
! !INTERFACE:
subroutine olprad
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Calculates the radial overlap integrals of the APW and local-orbital basis
!   functions. In other words, for atom $\alpha$, it computes integrals of the
!   form
!   $$ o^{\alpha}_{qp}=\int_0^{R_i}u^{\alpha}_{q;l_p}(r)v^{\alpha}_p(r)r^2dr $$
!   and
!   $$ o^{\alpha}_{pp'}=\int_0^{R_i}v^{\alpha}_p(r)v^{\alpha}_{p'}(r)r^2dr,
!    \quad l_p=l_{p'} $$
!   where $u^{\alpha}_{q;l}$ is the $q$th APW radial function for angular
!   momentum $l$; and $v^{\alpha}_p$ is the $p$th local-orbital radial function
!   and has angular momentum $l_p$.
!
! !REVISION HISTORY:
!   Created November 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ias,nr,nthd
integer ilo,jlo,l,io
! loop over atoms
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is,nr,ilo,jlo,l,io) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
! loop over local-orbitals
  do ilo=1,nlorb(is)
    l=lorbl(ilo,is)
!-------------------------------------!
!     APW-local-orbital integrals     !
!-------------------------------------!
    do io=1,apword(l,is)
      oalo(io,ilo,ias)=sum(apwfr(1:nr,1,io,l,ias)*lofr(1:nr,1,ilo,ias) &
       *wrmt(1:nr,is))
    end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
    do jlo=1,nlorb(is)
      if (lorbl(jlo,is) == l) then
        ololo(ilo,jlo,ias)=sum(lofr(1:nr,1,ilo,ias)*lofr(1:nr,1,jlo,ias) &
         *wrmt(1:nr,is))
      end if
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine
!EOC

