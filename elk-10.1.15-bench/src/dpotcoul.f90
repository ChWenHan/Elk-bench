
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dpotcoul
use modmain
use modphonon
implicit none
! local variables
integer np
! allocatable arrays
complex(8), allocatable :: gzfmt(:,:)
! solve the complex Poisson's equation in the muffin-tins
call genzvclmt(nrmt,nrmti,nrmtmax,rlmt,wprmt,npmtmax,drhomt,dvclmt)
! calculate the gradient of the nuclear potential
allocate(gzfmt(npmtmax,3))
call gradzvcln(isph,gzfmt)
! subtract gradient component corresponding to the phonon polarisation
np=npmt(isph)
dvclmt(1:np,iasph)=dvclmt(1:np,iasph)-gzfmt(1:np,ipph)
deallocate(gzfmt)
! solve Poisson's equation in the entire unit cell
call zpotcoul(0,nrmt,nrmti,npmt,nrmtmax,rlmt,ngridg,igfft,ngvec,gqc,gclgq, &
 ngvec,jlgqrmt,ylmgq,sfacgq,drhoir,npmtmax,dvclmt,dvclir)
end subroutine

