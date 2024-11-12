
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writemomru
use modmain
use modulr
implicit none
! local variables
integer is,ia,ias,ir
open(50,file='MOMENTRU.OUT',form='FORMATTED')
do ir=1,nqpt
  write(50,*)
  write(50,'("R-point (Cartesian coordinates) :")')
  write(50,'(3G18.10)') vrcu(:,ir)
  write(50,'("Moments :")')
  write(50,'(" interstitial",T30,": ",3G18.10)') momirru(1:ndmag,ir)
  write(50,'(" muffin-tins")')
  do is=1,nspecies
    write(50,'("  species : ",I4," (",A,")")') is,trim(spsymb(is))
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(50,'("   atom ",I4,T30,": ",3G18.10)') ia,mommtru(1:ndmag,ias,ir)
    end do
  end do
  write(50,'(" total moment",T30,": ",3G18.10)') momtotru(1:ndmag,ir)
end do
close(50)
end subroutine

