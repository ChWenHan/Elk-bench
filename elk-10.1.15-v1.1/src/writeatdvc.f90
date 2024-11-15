
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeatdvc
use modmain
implicit none
! local variables
integer is,ia
! write the atomic displacements and velocities in Cartesian coordinates to file
open(50,file='ATDVC.OUT',form='FORMATTED',action='WRITE')
do is=1,nspecies
  do ia=1,natoms(is)
    write(50,'(2I4,6G18.10)') is,ia,atdvc(:,:,ia,is)
  end do
end do
close(50)
end subroutine

