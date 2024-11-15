
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readepsw0
use modmain
use modphonon
implicit none
! local variables
logical exist
integer i,j
real(8) w
character(256) fname
do j=1,3
  do i=1,3
    write(fname,'("EPSILON_",2I1,".OUT")') i,j
    inquire(file=trim(fname),exist=exist)
    if (.not.exist) then
      write(*,*)
      write(*,'("Error(readepsw0): file not found :")')
      write(*,'(" ",A)') trim(fname)
      write(*,*)
      stop
    end if
    open(50,file=trim(fname),status='OLD',form='FORMATTED')
    read(50,*) w,epsw0(i,j)
    if (abs(w) > 1.d-8) then
      write(*,*)
      write(*,'("Error(readepsw0): first frequency should be zero")')
      write(*,'(" in file ",A)') trim(fname)
      write(*,*)
      stop
    end if
    close(50)
  end do
end do
! make dielectric tensor symmetric
do j=1,3
  do i=1,j-1
    epsw0(i,j)=0.5d0*(epsw0(i,j)+epsw0(j,i))
    epsw0(j,i)=epsw0(i,j)
  end do
end do
end subroutine

