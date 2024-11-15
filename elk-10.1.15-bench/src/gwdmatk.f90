
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gwdmatk(ik)
use modmain
use modgw
use modomp
implicit none
! arguments
integer, intent(in) :: ik
! local variables
integer ist,jst,iw,i,nthd
real(8) e,t1
complex(8) z1
! automatic arrays
complex(8) gs(nstsv),g(nstsv,nstsv),ge(4,nstsv,nstsv)
complex(8) evecsv(nstsv,nstsv),d(nstsv,nstsv),a(nstsv,nstsv)
! allocatable arrays
complex(8), allocatable :: se(:,:,:)
! external functions
complex(8), external :: gwtails
! read the self-energy from file
allocate(se(nstsv,nstsv,0:nwfm))
call getgwsefm(ik,se)
! zero the density matrix
d(1:nstsv,1:nstsv)=0.d0
! loop over fermionic Matsubara frequencies
call holdthd(nwfm+1,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(gs,g,ist,e,z1,i) &
!$OMP REDUCTION(+:d) &
!$OMP NUM_THREADS(nthd)
do iw=0,nwfm
! compute the diagonal matrix G_s
  do ist=1,nstsv
    e=evalsv(ist,ik)-efermi
    gs(ist)=1.d0/(wfm(iw)-e)
  end do
! compute 1 - G_s Σ
  do ist=1,nstsv
    z1=-gs(ist)
    g(ist,1:nstsv)=z1*se(ist,1:nstsv,iw)
    g(ist,ist)=g(ist,ist)+1.d0
  end do
! invert this matrix
  call zminv(nstsv,g)
! compute G = (1 - G_s Σ)⁻¹ G_s
  do jst=1,nstsv
    z1=gs(jst)
    g(1:nstsv,jst)=g(1:nstsv,jst)*z1
  end do
! add to the density matrix
  d(1:nstsv,1:nstsv)=d(1:nstsv,1:nstsv)+g(1:nstsv,1:nstsv)
! store the Green's function at the end point frequencies
  i=0
  if (iw == 0) i=1
  if (iw == 1) i=2
  if (iw == nwfm-1) i=3
  if (iw == nwfm) i=4
  if (i /= 0) ge(i,1:nstsv,1:nstsv)=g(1:nstsv,1:nstsv)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! add the Matsubara tails analytically
do jst=1,nstsv
  do ist=1,nstsv
    d(ist,jst)=d(ist,jst)+gwtails(ge(:,ist,jst))
  end do
end do
! multiply by 1/β
t1=kboltz*tempk
d(1:nstsv,1:nstsv)=t1*d(1:nstsv,1:nstsv)
deallocate(se)
! make density matrix Hermitian
do ist=1,nstsv
  do jst=1,ist-1
    z1=0.5d0*(d(ist,jst)+conjg(d(jst,ist)))
    d(ist,jst)=z1
    d(jst,ist)=conjg(z1)
  end do
  d(ist,ist)=dble(d(ist,ist))
end do
! diagonalise the density matrix for the natural orbitals and occupation numbers
call eveqnzh(nstsv,nstsv,d,occsv(:,ik))
occsv(1:nstsv,ik)=occsv(1:nstsv,ik)*occmax
! get the second-variational eigenvectors from file
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! apply unitary transformation to the third-variational states so that they
! refer to the first-variational basis
call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,d,nstsv,zzero,a,nstsv)
! write the density matrix to file as second-variational eigenvectors
call putevecsv(filext,ik,a)
end subroutine

