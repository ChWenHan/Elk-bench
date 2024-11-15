
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradzf(zfmt,zfir,gzfmt,gzfir)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(in) :: zfmt(npmtmax,natmtot),zfir(ngtot)
complex(8), intent(out) :: gzfmt(npmtmax,natmtot,3),gzfir(ngtot,3)
! local variables
integer is,ias,ld,i
integer ig,ifg,nthd
complex(8) z1
! allocatable arrays
complex(8), allocatable :: zfft(:)
! muffin-tin gradient
ld=npmtmax*natmtot
call holdthd(natmtot+1,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfft,is,i,ig,ifg,z1) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call gradzfmt(nrmt(is),nrmti(is),rlmt(:,-1,is),wcrmt(:,:,is),zfmt(:,ias),ld, &
   gzfmt(1,ias,1))
end do
!$OMP END DO NOWAIT
! interstitial gradient
!$OMP SINGLE
allocate(zfft(ngtot))
zfft(1:ngtot)=zfir(1:ngtot)
call zfftifc(3,ngridg,-1,zfft)
do i=1,3
  do ig=1,ngvec
    ifg=igfft(ig)
    z1=zfft(ifg)
    gzfir(ifg,i)=vgc(i,ig)*cmplx(-aimag(z1),dble(z1),8)
  end do
  gzfir(igfft(ngvec+1:ngtot),i)=0.d0
  call zfftifc(3,ngridg,1,gzfir(:,i))
end do
deallocate(zfft)
!$OMP END SINGLE
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

