
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine curlrvf(rvfmt,rvfir,curlmt,curlir)
use modmain
use modomp
implicit none
! arguments
real(8), intent(in) :: rvfmt(npmtmax,natmtot,3),rvfir(ngtot,3)
real(8), intent(out) :: curlmt(npmtmax,natmtot,3),curlir(ngtot,3)
! local variables
integer is,ias,np,i,nthd
! allocatable arrays
real(8), allocatable :: grfmt(:,:,:,:),grfir(:,:,:)
allocate(grfmt(npmtmax,natmtot,3,3),grfir(ngtot,3,3))
! compute the gradients
call holdthd(3,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,ias,is,np) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do i=1,3
  call gradrf(rvfmt(:,:,i),rvfir(:,i),grfmt(:,:,:,i),grfir(:,:,i))
end do
!$OMP END DO
! determine the muffin-tin and interstitial curl
!$OMP SECTIONS
!$OMP SECTION
do ias=1,natmtot
  is=idxis(ias)
  np=npmt(is)
  curlmt(1:np,ias,1)=grfmt(1:np,ias,2,3)-grfmt(1:np,ias,3,2)
  curlmt(1:np,ias,2)=grfmt(1:np,ias,3,1)-grfmt(1:np,ias,1,3)
  curlmt(1:np,ias,3)=grfmt(1:np,ias,1,2)-grfmt(1:np,ias,2,1)
end do
!$OMP SECTION
curlir(1:ngtot,1)=grfir(1:ngtot,2,3)-grfir(1:ngtot,3,2)
curlir(1:ngtot,2)=grfir(1:ngtot,3,1)-grfir(1:ngtot,1,3)
curlir(1:ngtot,3)=grfir(1:ngtot,1,2)-grfir(1:ngtot,2,1)
!$OMP END SECTIONS
!$OMP END PARALLEL
call freethd(nthd)
deallocate(grfmt,grfir)
end subroutine

