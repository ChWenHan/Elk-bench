
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetdforces
use modmain
use modtddft
implicit none
! local variables
integer is,ia,ias
! write the total force on each atom
open(50,file='FORCETOT_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(I8,G18.10)') itimes,times(itimes)
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,'(2I4,3G18.10)') is,ia,forcetot(1:3,ias)
  end do
end do
close(50)
! write the maximum force magnitude over all atoms
open(50,file='FORCEMAX_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(2G18.10)') times(itimes),forcemax
close(50)
end subroutine

