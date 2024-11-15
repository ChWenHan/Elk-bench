
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine batchdv
use modmain
use moddftu
use modvars
use modmpi
implicit none
! no increment for first task
if (itask <= 1) return
! only increment for ground-state tasks
if (all(task /= [0,1,2,3])) return
! increment selected variables
if (any(dngridk(:) /= 0)) then
  ngridk(:)=ngridk(:)+dngridk(:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented ngridk")')
end if
if (dlmaxapw /= 0) then
  lmaxapw=lmaxapw+dlmaxapw
  if (mp_mpi) write(*,'("Info(batchdv): incremented lmaxapw")')
end if
if (dlmaxo /= 0) then
  lmaxo=lmaxo+dlmaxo
  if (mp_mpi) write(*,'("Info(batchdv): incremented lmaxo")')
end if
if (any(davec(:,:) /= 0.d0)) then
  avec(:,:)=avec(:,:)+davec(:,:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented avec")')
end if
if (any(datposl(:,:,:) /= 0.d0)) then
  atposl(:,:,:)=atposl(:,:,:)+datposl(:,:,:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented atposl")')
end if
if (drgkmax /= 0.d0) then
  rgkmax=rgkmax+drgkmax
  if (mp_mpi) write(*,'("Info(batchdv): incremented rgkmax")')
end if
if (dgmaxvr /= 0.d0) then
  gmaxvr=gmaxvr+dgmaxvr
  if (mp_mpi) write(*,'("Info(batchdv): incremented gmaxvr")')
end if
if (dnempty0 /= 0.d0) then
  nempty0=nempty0+dnempty0
  if (mp_mpi) write(*,'("Info(batchdv): incremented nempty")')
end if
if (dchgexs /= 0.d0) then
  chgexs=chgexs+dchgexs
  if (mp_mpi) write(*,'("Info(batchdv): incremented chgexs")')
end if
if (dsxcscf /= 0.d0) then
  sxcscf=sxcscf+dsxcscf
  if (mp_mpi) write(*,'("Info(batchdv): incremented sxcscf")')
end if
if (dnrmtscf /= 0.d0) then
  nrmtscf=nrmtscf+dnrmtscf
  if (mp_mpi) write(*,'("Info(batchdv): incremented nrmtscf")')
end if
if (any(dudufix(:) /= 0.d0)) then
  udufix(:)=udufix(:)+dudufix(:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented udufix")')
end if
if (any(dmomfix(:) /= 0.d0)) then
  momfix(:)=momfix(:)+dmomfix(:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented momfix")')
end if
if (any(dbfieldc0(:) /= 0.d0)) then
  bfieldc0(:)=bfieldc0(:)+dbfieldc0(:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented dbfieldc0")')
end if
if (any(dvqlss(:) /= 0.d0)) then
  vqlss(:)=vqlss(:)+dvqlss(:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented vqlss")')
end if
if (any(dafieldc(:) /= 0.d0)) then
  afieldc(:)=afieldc(:)+dafieldc(:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented afieldc")')
end if
if (any(dafspc(:,:) /= 0.d0)) then
  afspc(:,:)=afspc(:,:)+dafspc(:,:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented afspc")')
end if
end subroutine

