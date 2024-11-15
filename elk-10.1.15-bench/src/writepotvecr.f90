
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writepotvecr(fname,potvecr)
use modmain
use modulr
implicit none
!
character(*), intent(in) :: fname
real(8), intent(out) :: potvecr(3,nqpt)
! local variables
integer i1,i2,i3,ir
! write the real-space potential to file
open(50,file=trim(fname),form='FORMATTED')
write(50,'(3I6," : grid size")') ngridq
ir=0
do i3=1,ngridq(3)
  do i2=1,ngridq(2)
    do i1=1,ngridq(1)
      ir=ir+1
      write(50,'(3I6,3G18.10)') i1,i2,i3,potvecr(1:3,ir)
    end do
  end do
end do
close(50)
end subroutine

