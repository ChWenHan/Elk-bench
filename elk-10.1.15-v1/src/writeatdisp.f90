
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeatdisp
use modmain
use modtddft
implicit none
! local variables
integer is,ia
real(8) vl(3)
open(50,file='ATDISPL_TD.OUT',form='FORMATTED',position='APPEND')
open(51,file='ATDISPC_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(I8,G18.10)') itimes,times(itimes)
write(51,'(I8,G18.10)') itimes,times(itimes)
do is=1,nspecies
  do ia=1,natoms(is)
    call r3mv(ainv,atdvc(:,0,ia,is),vl)
    write(50,'(2I4,3G18.10)') is,ia,vl(:)
    write(51,'(2I4,3G18.10)') is,ia,atdvc(:,0,ia,is)
  end do
end do
close(50)
close(51)
end subroutine

