
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initph
use modmain
use modphonon
implicit none
! allocate global arrays
if (allocated(dynq)) deallocate(dynq)
allocate(dynq(nbph,nbph,nqpt))
if (allocated(dynr)) deallocate(dynr)
allocate(dynr(nbph,nbph,nqptnr))
if (allocated(wphq)) deallocate(wphq)
allocate(wphq(nbph,nqpt))
! read in the dynamical matrices
call readdyn
! apply the acoustic sum rule
call sumrule
! if the non-analytic term is required then read in the Born effective charge
! tensor for each atom as well as the static dielectric tensor
if (allocated(bec)) deallocate(bec)
if (tphnat) then
  allocate(bec(3,3,natmtot))
  call readbec
  call readepsw0
end if
! Fourier transform the dynamical matrices to real-space
call dynqtor(dynq,dynr)
end subroutine

