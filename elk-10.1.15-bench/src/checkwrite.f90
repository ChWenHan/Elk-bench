
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine checkwrite(twrite)
use modmain
use moddelf
implicit none
! arguments
logical, intent(out) :: twrite
! check for WRITE file
inquire(file='WRITE',exist=twrite)
if (twrite) then
  write(*,'("Info(checkwrite): WRITE file exists")')
! delete the WRITE file
  call delfile('WRITE')
end if
end subroutine

