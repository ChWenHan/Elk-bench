
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readbec
use modmain
use modphonon
implicit none
! local variables
logical exist
integer is,ia,ias,ip,jp
character(256) fext
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ip=1,3
      call becfext(is,ia,ip,fext)
      inquire(file='BEC'//trim(fext),exist=exist)
      if (.not.exist) then
        write(*,*)
        write(*,'("Error(readbec): file not found :")')
        write(*,'(A)') ' BEC'//trim(fext)
        write(*,*)
        stop
      end if
      open(50,file='BEC'//trim(fext),status='OLD',form='FORMATTED')
      do jp=1,3
        read(50,*) bec(ip,jp,ias)
      end do
      close(50)
    end do
  end do
end do
end subroutine

