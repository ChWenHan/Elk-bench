
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genidxthc
use modmain
use modtdhfc
implicit none
! local variables
integer ik,nst,nmax,ist
! determine the maximum number of states within energy window over all k-points
nmax=0
do ik=1,nkpt
  nst=0
  do ist=1,nstsv
    if (abs(evalsv(ist,ik)-efermi) < ecutthc) nst=nst+1
  end do
  nmax=max(nmax,nst)
end do
if (nmax == 0) then
  write(*,*)
  write(*,'("Error(genidxthc): no states within energy window ecutthc")')
  write(*,*)
  stop
end if
! allocate global arrays
if (allocated(nstthc)) deallocate(nstthc)
allocate(nstthc(nkpt))
if (allocated(idxthc)) deallocate(idxthc)
allocate(idxthc(nmax,nkpt))
! determine the number of and index to used states
do ik=1,nkpt
  nst=0
  do ist=1,nstsv
    if (abs(evalsv(ist,ik)-efermi) < ecutthc) then
      nst=nst+1
      idxthc(nst,ik)=ist
    end if
  end do
  nstthc(ik)=nst
end do
end subroutine

