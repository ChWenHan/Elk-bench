
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine energynn
use modmain
implicit none
! local variables
integer is,ias,nr,nri
integer iro,i0,i1
real(8) t1
! allocatable arrays
complex(8), allocatable :: zvclmt(:,:),zvclir(:),zrhoir(:)
allocate(zvclmt(npmtmax,natmtot))
! generate the nuclear monopole potentials
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  iro=nri+1
  i1=lmmaxi*(nri-1)+1
  zvclmt(1:npmt(is),ias)=0.d0
  zvclmt(1:i1:lmmaxi,ias)=vcln(1:nri,is)
  i0=i1+lmmaxi
  i1=lmmaxo*(nr-iro)+i0
  zvclmt(i0:i1:lmmaxo,ias)=vcln(iro:nr,is)
end do
allocate(zrhoir(ngtot),zvclir(ngtot))
! set the interstitial density to zero
zrhoir(1:ngtot)=0.d0
! solve the complex Poisson's equation
call zpotcoul(0,nrmt,nrmti,npmt,nrmtmax,rlmt,ngridg,igfft,ngvec,gc,gclg,ngvec, &
 jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)
! compute the nuclear-nuclear energy
engynn=0.d0
do ias=1,natmtot
  is=idxis(ias)
  t1=(dble(zvclmt(1,ias))-vcln(1,is))*y00
  engynn=engynn+spzn(is)*t1
end do
engynn=0.5d0*engynn
deallocate(zvclmt,zvclir,zrhoir)
end subroutine

