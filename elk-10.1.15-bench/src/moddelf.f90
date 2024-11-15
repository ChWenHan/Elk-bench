
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module moddelf

use modphonon
use modramdisk
use modmpi

contains

subroutine delfile(fname)
implicit none
! arguments
character(*), intent(in) :: fname
! local variables
integer ios
open(40,file=fname,iostat=ios)
close(40,status='DELETE',iostat=ios)
end subroutine

subroutine delfiles(evec,devec,eval,occ,pmat,epsi)
implicit none
! arguments
logical, optional, intent(in) :: evec,devec,eval,occ,pmat,epsi
! local variables
character(256) fext,fname
if (present(evec)) then
! delete the first-variational eigenvector file
  fname=trim(scrpath)//'EVECFV'//trim(filext)
  if (mp_mpi) call delfile(fname)
  if (ramdisk) call delfrd(fname)
! delete the second-variational eigenvector file
  fname=trim(scrpath)//'EVECSV'//trim(filext)
  if (mp_mpi) call delfile(fname)
  if (ramdisk) call delfrd(fname)
end if
if (present(devec)) then
! construct the dynamical matrix file extension
  call dynfext(iqph,isph,iaph,ipph,fext)
! delete the eigenvector derivative files
  fname=trim(scrpath)//'DEVECFV'//trim(fext)
  if (mp_mpi) call delfile(fname)
  if (ramdisk) call delfrd(fname)
  fname=trim(scrpath)//'DEVECSV'//trim(fext)
  if (mp_mpi) call delfile(fname)
  if (ramdisk) call delfrd(fname)
end if
if (present(eval)) then
! delete the first-variational eigenvalue file
  fname='EVALFV'//trim(filext)
  if (mp_mpi) call delfile(fname)
  if (ramdisk) call delfrd(fname)
! delete the second-variational eigenvalue file
  if (mp_mpi) call delfile('EVALSV'//trim(filext))
end if
if (present(occ)) then
! delete the occupation number file
  if (mp_mpi) call delfile('OCCSV'//trim(filext))
end if
if (present(pmat)) then
! delete the momentum matrix elements file
  if (mp_mpi) call delfile('PMAT.OUT')
  if (ramdisk) call delfrd('PMAT.OUT')
end if
if (present(epsi)) then
! delete the inverse dielectric function file
  if (mp_mpi) call delfile('EPSINV.OUT')
end if
end subroutine

end module

