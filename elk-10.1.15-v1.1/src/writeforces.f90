
! Copyright (C) 2004-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeforces(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer is,ia,ias
write(fnum,*)
write(fnum,'("Forces :")')
do is=1,nspecies
  write(fnum,'(" species : ",I4," (",A,")")') is,trim(spsymb(is))
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(fnum,'("  atom : ",I4)') ia
    write(fnum,'("   Hellmann-Feynman",T30,": ",3F14.8)') forcehf(:,ias)
    write(fnum,'("   IBS",T30,": ",3F14.8)') forcetot(:,ias)-forcehf(:,ias)
    write(fnum,'("   total force",T30,": ",3F14.8)') forcetot(:,ias)
    write(fnum,'("   total magnitude",T30,": ",F14.8)') norm2(forcetot(:,ias))
  end do
end do
end subroutine

