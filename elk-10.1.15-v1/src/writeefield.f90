
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeefield(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer is,ia,ias,i
real(8) efc(3),t1
if (natmtot == 0) return
! determine the average electric field in each muffin-tin
call efieldmt
! write the electric fields to file
write(fnum,*)
write(fnum,'("Average electric field in each muffin-tin")')
do is=1,nspecies
  write(fnum,'("  species : ",I4," (",A,")")') is,trim(spsymb(is))
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(fnum,'("   atom ",I4,T30,": ",3G18.10)') ia,efcmt(:,ias)
  end do
end do
! compute the average electric field
do i=1,3
  efc(i)=sum(efcmt(i,1:natmtot))/dble(natmtot)
end do
write(fnum,*)
write(fnum,'("Average of muffin-tin electric fields :")')
write(fnum,'(3G18.10)') efc
t1=norm2(efc(:))
write(fnum,'(" magnitude : ",G18.10)') t1
write(fnum,'("  volts/nanometer : ",G18.10)') t1*ef_si/1.d9
end subroutine

