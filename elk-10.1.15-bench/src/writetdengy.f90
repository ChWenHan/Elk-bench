
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetdengy
use modmain
use modtddft
implicit none
open(50,file='TOTENERGY_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(2G18.10)') times(itimes),engytot
close(50)
end subroutine

