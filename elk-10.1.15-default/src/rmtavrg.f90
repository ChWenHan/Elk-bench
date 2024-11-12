
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rmtavrg
use modmain
implicit none
! local variables
integer is,i
real(8) ra
if (nspecies <= 1) return
do i=1,mrmtav
! average muffin-tin radius
  ra=sum(rmt(1:nspecies))/nspecies
! replace each muffin-tin radius with half itself plus the average
  do is=1,nspecies
    rmt(is)=0.5d0*(rmt(is)+ra)
  end do
end do
end subroutine

