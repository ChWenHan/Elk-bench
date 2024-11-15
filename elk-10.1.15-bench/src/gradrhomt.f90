
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradrhomt
use modmain
use modphonon
implicit none
! local variables
integer nr,nri,np
! allocatable arrays
complex(8), allocatable :: zfmt(:),grhomt(:,:)
allocate(zfmt(npmtmax),grhomt(npmtmax,3))
! add gradient contribution from rigid shift of muffin-tin
nr=nrmt(isph)
nri=nrmti(isph)
np=npmt(isph)
! convert the density to complex spherical harmonic expansion
call rtozfmt(nr,nri,rhomt(:,iasph),zfmt)
! compute the gradient
call gradzfmt(nr,nri,rlmt(:,-1,isph),wcrmt(:,:,isph),zfmt,npmtmax,grhomt)
! store density derivative for displaced atom
if (tlast) drhomt(1:np,natmtot+1)=drhomt(1:np,iasph)
! subtract from the density derivative
drhomt(1:np,iasph)=drhomt(1:np,iasph)-grhomt(1:np,ipph)
deallocate(zfmt,grhomt)
end subroutine

