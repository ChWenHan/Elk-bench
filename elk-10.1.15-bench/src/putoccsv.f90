
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putoccsv(fext,ik,occsvp)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
real(8), intent(in) :: occsvp(nstsv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(1:3,ik),nstsv,occsvp
!$OMP CRITICAL(u208)
open(208,file='OCCSV'//trim(fext),form='UNFORMATTED',access='DIRECT',recl=recl)
write(208,rec=ik) vkl(1:3,ik),nstsv,occsvp
close(208)
!$OMP END CRITICAL(u208)
end subroutine

