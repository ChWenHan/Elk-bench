
! Copyright (C) 2014 D. Ernsting, S. Dugdale and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine reademd(emds)
use modmain
use modpw
use modtest
implicit none
! arguments
real(4), intent(out) :: emds(nhkmax,nkpt)
! local variables
integer ik,recl,nhk_
real(8) vkl_(3),t1
! allocatable arrays
real(8), allocatable :: emd(:)
allocate(emd(nhkmax))
! find the record length
inquire(iolength=recl) vkl_,nhk_,emd
open(250,file='EMD.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
do ik=1,nkpt
  read(250,rec=ik) vkl_,nhk_,emd
  t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
  if (t1 > epslat) then
    write(*,*)
    write(*,'("Error(reademd): differing vectors for k-point ",I8)') ik
    write(*,'(" current : ",3G18.10)') vkl(:,ik)
    write(*,'(" EMD.OUT : ",3G18.10)') vkl_
    write(*,*)
    stop
  end if
  if (nhk(1,ik) /= nhk_) then
    write(*,*)
    write(*,'("Error(reademd): differing nhk for k-point ",I8)') ik
    write(*,'(" current : ",I8)') nhk(1,ik)
    write(*,'(" EMD.OUT : ",I8)') nhk_
    write(*,*)
    stop
  end if
! store momentum density in single-precision array
  emds(:,ik)=real(emd(:))
end do
close(250)
! write the EMD of the last k-point to test file if required
call writetest(170,'electron momentum density',nv=nhkmax,tol=1.d-4,rva=emd)
deallocate(emd)
end subroutine

