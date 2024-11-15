
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initoep
use modmain
implicit none
! local variables
integer is,ist,nc
! find maximum core states over all species
ncrmax=0
do is=1,nspecies
  nc=0
  do ist=1,nstsp(is)
    if (spcore(ist,is)) nc=nc+2*ksp(ist,is)
  end do
  ncrmax=max(ncrmax,nc)
end do
! allocate the exchange potential and magnetic field
if (allocated(vxmt)) deallocate(vxmt)
allocate(vxmt(npcmtmax,natmtot))
if (allocated(vxir)) deallocate(vxir)
allocate(vxir(ngtot))
if (spinpol) then
  if (allocated(bxmt)) deallocate(bxmt)
  allocate(bxmt(npcmtmax,natmtot,ndmag))
  if (allocated(bxir)) deallocate(bxir)
  allocate(bxir(ngtot,ndmag))
end if
! allocate the OEP residual functions
allocate(dvxmt(npcmtmax,natmtot),dvxir(ngtot))
if (spinpol) then
  allocate(dbxmt(npcmtmax,natmtot,ndmag),dbxir(ngtot,ndmag))
end if
! set initial step size for iterative method
tauoep=tau0oep
end subroutine

