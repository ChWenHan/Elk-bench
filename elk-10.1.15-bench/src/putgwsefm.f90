
! Copyright (C) 2017 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine putgwsefm(ik,se)
use modmain
use modgw
use modramdisk
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: se(nstsv,nstsv,0:nwfm)
! local variables
integer recl
!$OMP CRITICAL(u280)
! write to RAM disk if required
if (ramdisk) then
  call putrd('GWSEFM.OUT',ik,v1=vkl(1:3,ik),n1=nstsv,n2=nwfm, &
   nzv=nstsv*nstsv*(nwfm+1),zva=se)
end if
! write to disk if required
if (wrtdsk) then
! find the record length
  inquire(iolength=recl) vkl(1:3,ik),nstsv,nwfm,se
  open(280,file='GWSEFM.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
  write(280,rec=ik) vkl(1:3,ik),nstsv,nwfm,se
  close(280)
end if
!$OMP END CRITICAL(u280)
end subroutine

