
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writebox(fnum,str)
implicit none
! arguments
integer, intent(in) :: fnum
character(*), intent(in) :: str
! local variables
integer l
l=len(str)+2
write(fnum,*)
write(fnum,'("┌",A,"┐")') repeat("─",l)
write(fnum,'("│ ",A," │")') str
write(fnum,'("└",A,"┘")') repeat("─",l)
end subroutine

