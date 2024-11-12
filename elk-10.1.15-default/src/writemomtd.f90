
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writemomtd
use modmain
use modtddft
implicit none
! local variables
integer is,ia,ias
! write the total spin moment
open(50,file='MOMENT_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(4G18.10)') times(itimes),momtot(1:ndmag)
close(50)
! write the total spin moment magnitude
open(50,file='MOMENTM_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(2G18.10)') times(itimes),momtotm
close(50)
! write muffin-tin moments
open(50,file='MOMENTMT_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(G18.10)') times(itimes)
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,'(2I4,3G18.10)') is,ia,mommt(1:ndmag,ias)
  end do
end do
write(50,*)
close(50)
! write interstitial moment
open(50,file='MOMENTIR_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(4G18.10)') times(itimes),momir(1:ndmag)
close(50)
end subroutine

