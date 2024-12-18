
! Copyright (C) 2015 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tdtemp(occsvp)
use modmain
use modtddft
implicit none
! arguments
real(8), intent(in) :: occsvp(nstsv,nkpt)
! local variables
integer, parameter :: maxit=1000
integer ik,ist,it
real(8), parameter :: eps=1.d-6
real(8) sw,dsw,sm,sp
real(8) x,t1,t2,t3
! external functions
real(8), external :: sdelta_fd,stheta_fd
! initial smearing width
sw=1.d-6
! initial smearing width step size
dsw=1.d-6
sp=0.d0
do it=1,maxit
  t1=1.d0/sw
  sm=0.d0
  do ik=1,nkpt
    do ist=1,nstsv
      x=(efermi-evalsv(ist,ik))*t1
      t2=occmax*stheta_fd(x)
      t3=sdelta_fd(x)*x*t1
      sm=sm+(occsvp(ist,ik)-t2)*t3
    end do
  end do
  if ((sm*sp) < 0.d0) then
    dsw=-0.5d0*dsw
    if (abs(dsw) < eps) goto 10
  else
    dsw=1.5d0*dsw
  end if
  sp=sm
  sw=sw+dsw
  if ((sw < 0.d0).or.(sw > 1.d6)) exit
end do
write(*,*)
write(*,'("Warning(tdtemp): could not estimate effective temperature")')
return
10 continue
! write effective temperature to file
t1=sw/kboltz
open(50,file='TDTEMP.OUT',form='FORMATTED',position='APPEND')
write(50,'(2G18.10)') times(itimes),t1
close(50)
end subroutine

