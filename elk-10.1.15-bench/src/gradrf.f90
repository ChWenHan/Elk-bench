
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradrf(rfmt,rfir,grfmt,grfir)
use modmain
use modomp
implicit none
! arguments
real(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
real(8), intent(out) :: grfmt(npmtmax,natmtot,3),grfir(ngtot,3)
! local variables
integer is,ias,ld,i
integer ig,ifg,nthd
complex(8) z1
! allocatable arrays
complex(8), allocatable :: zfft1(:),zfft2(:)
! muffin-tin gradient
ld=npmtmax*natmtot
call holdthd(natmtot+1,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfft1,zfft2) &
!$OMP PRIVATE(is,i,ifg,ig,z1) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call gradrfmt(nrmt(is),nrmti(is),rlmt(:,-1,is),wcrmt(:,:,is),rfmt(:,ias),ld, &
   grfmt(1,ias,1))
end do
!$OMP END DO NOWAIT
! interstitial gradient
!$OMP SINGLE
allocate(zfft1(nfgrz),zfft2(nfgrz))
call rzfftifc(3,ngridg,-1,rfir,zfft1)
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    if (ig <= ngvec) then
      z1=zfft1(ifg)
      zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(z1),dble(z1),8)
    else
      zfft2(ifg)=0.d0
    end if
  end do
  call rzfftifc(3,ngridg,1,grfir(:,i),zfft2)
end do
deallocate(zfft1,zfft2)
!$OMP END SINGLE
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

