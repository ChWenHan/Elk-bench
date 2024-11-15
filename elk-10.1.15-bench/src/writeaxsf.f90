
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeaxsf
use modmain
use modtddft
implicit none
! local variables
integer is,ia
! write header at first time step
if (itimes <= 1) then
  open(50,file='crystal.axsf',form='FORMATTED')
  write(50,'("ANIMSTEPS ",I8)') (ntimes-2)/ntsforce+1
  write(50,'("CRYSTAL")')
  write(50,'("PRIMVEC")')
  write(50,'(3G18.10)') avec(:,1)*br_ang
  write(50,'(3G18.10)') avec(:,2)*br_ang
  write(50,'(3G18.10)') avec(:,3)*br_ang
  close(50)
end if
open(50,file='crystal.axsf',form='FORMATTED',position='APPEND')
write(50,*)
write(50,'("PRIMCOORD ",I8)') (itimes-1)/ntsforce+1
write(50,'(2I8)') natmtot,1
do is=1,nspecies
  do ia=1,natoms(is)
    write(50,'(A,3G18.10)') trim(spsymb(is)),atposc(:,ia,is)*br_ang
  end do
end do
close(50)
end subroutine

