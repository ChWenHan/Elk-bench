
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readpotvecr(fname,potvecr)
use modmain
use modulr
implicit none
!
character(*), intent(in) :: fname
real(8), intent(out) :: potvecr(3,nqpt)
! local variables
integer i1,i2,i3,ir
integer ngridq_(3),i1_,i2_,i3_
! read the real-space external vector potential from file
open(50,file=trim(fname),form='FORMATTED')
read(50,*) ngridq_(:)
if (any(ngridq(:) /= ngridq_(:))) then
  write(*,*)
  write(*,'("Error(readpotvecr): differing ngridq")')
  write(*,'(" current  : ",3I6)') ngridq
  write(*,'(A," : ",3I6)') fname, ngridq_
  write(*,*)
  stop
end if
ir=0
do i3=1,ngridq(3)
  do i2=1,ngridq(2)
    do i1=1,ngridq(1)
      ir=ir+1
      read(50,*) i1_,i2_,i3_,potvecr(1,ir),potvecr(2,ir),potvecr(3,ir)
      if ((i1 /= i1_).or.(i2 /= i2_).or.(i3 /= i3_)) then
        write(*,*)
        write(*,'("Error(readpotvecr): differing i1, i2 or i3")')
        write(*,'(" current  : ",3I6)') i1,i2,i3
        write(*,'(A," : ",3I6)') fname, i1_,i2_,i3_
        write(*,*)
        stop
      end if
    end do
  end do
end do
close(50)
end subroutine

