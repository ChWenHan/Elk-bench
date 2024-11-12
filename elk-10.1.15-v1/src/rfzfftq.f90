
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfzfftq(sgn,nf,ngt,rfmt,rfir,zfmt,zfir)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: sgn,nf,ngt
real(8), intent(inout) :: rfmt(npcmtmax,natmtot,nf,nqpt)
real(8), intent(inout) :: rfir(ngt,nf,nqpt)
complex(8), intent(inout) :: zfmt(npcmtmax,natmtot,nf,nfqrz)
complex(8), intent(inout) :: zfir(ngt,nf,nfqrz)
! local variables
integer jf,is,ias,npc,i,nthd
! automatic arrays
real(8) r(nqpt)
complex(8) z(nfqrz)
call holdthd(npcmtmax+ngt,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(r,z,jf,ias,is,npc,i) &
!$OMP NUM_THREADS(nthd)
if (sgn == -1) then
! loop over the number of functions
  do jf=1,nf
! Fourier transform the muffin-tin function
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
!$OMP DO
      do i=1,npc
        r(1:nqpt)=rfmt(i,ias,jf,1:nqpt)
        call rzfftifc(3,ngridq,-1,r,z)
        zfmt(i,ias,jf,1:nfqrz)=z(1:nfqrz)
      end do
!$OMP END DO NOWAIT
    end do
! Fourier transform the interstitial function
!$OMP DO
    do i=1,ngt
      r(1:nqpt)=rfir(i,jf,1:nqpt)
      call rzfftifc(3,ngridq,-1,r,z)
      zfir(i,jf,1:nfqrz)=z(1:nfqrz)
    end do
!$OMP END DO NOWAIT
! end loop over number of functions
  end do
else
! loop over the number of functions
  do jf=1,nf
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
!$OMP DO
      do i=1,npc
        z(1:nfqrz)=zfmt(i,ias,jf,1:nfqrz)
        call rzfftifc(3,ngridq,1,r,z)
        rfmt(i,ias,jf,1:nqpt)=r(1:nqpt)
      end do
!$OMP END DO NOWAIT
    end do
!$OMP DO
    do i=1,ngt
      z(1:nfqrz)=zfir(i,jf,1:nfqrz)
      call rzfftifc(3,ngridq,1,r,z)
      rfir(i,jf,1:nqpt)=r(1:nqpt)
    end do
!$OMP END DO NOWAIT
! end loop over number of functions
  end do
end if
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

