
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetddft
use modmain
use modtddft
use moddftu
use moddelf
implicit none
! local variables
integer is,ia,ias
character(32) fext
! allocatable arrays
real(8), allocatable :: rvfmt(:,:,:),rvfir(:,:)
! file extension
write(fext,'("_TS",I8.8,".OUT")') itimes
! delete files at first time step
if (itimes <= 1) then
  call delfile('TOTENERGY_TD.OUT')
  call delfile('CHARGEMT_TD.OUT')
  call delfile('CHARGEIR_TD.OUT')
  if (spinpol) then
    call delfile('MOMENT_TD.OUT')
    call delfile('MOMENTM_TD.OUT')
    call delfile('MOMENTMT_TD.OUT')
    call delfile('MOMENTIR_TD.OUT')
  end if
  call delfile('JTOT_TD.OUT')
  call delfile('JTOTM_TD.OUT')
  if (tddos) call delfile('TDTEMP.OUT')
  if (tdlsj) call delfile('TDLSJ.OUT')
  if (tforce) call delfile('FORCETOT_TD.OUT')
  if (tatdisp) then
    call delfile('ATDISPL_TD.OUT')
    call delfile('ATDISPC_TD.OUT')
  end if
  if (tafindt) call delfile('AFIND_TD.OUT')
end if
! total energy
call writetdengy
! muffin-tin charges
open(50,file='CHARGEMT_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(G18.10)') times(itimes)
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,'(2I4,G18.10)') is,ia,chgmt(ias)
  end do
end do
write(50,*)
close(50)
! interstitial charge
open(50,file='CHARGEIR_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(2G18.10)') times(itimes),chgir
close(50)
! write spin moments to file if required
if (spinpol) call writemomtd
! total current
open(50,file='JTOT_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(4G18.10)') times(itimes),jtot(:)
close(50)
! total current magnitude
open(50,file='JTOTM_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(4G18.10)') times(itimes),jtotm
close(50)
! write the time-dependent atomic forces
if (tforce.and.ttsforce) then
  call writetdforces
  open(50,file='FORCES'//trim(fext),form='FORMATTED',action='WRITE')
  call writeforces(50)
  close(50)
end if
! write the time-dependent atomic displacements
if (tatdisp.and.ttsforce) call writeatdisp
! write the time-dependent induced A-field
if (tafindt) then
  open(50,file='AFIND_TD.OUT',form='FORMATTED',position='APPEND')
  write(50,'(4G18.10)') times(itimes),afindt(:,0)
  close(50)
end if
! write optional quantities
if (.not.ttswrite) return
! charge density in 1D
if (tdrho1d) then
  open(50,file='RHO1D'//trim(fext),form='FORMATTED',action='WRITE')
  open(51,file='RHOLINES.OUT',form='FORMATTED',action='WRITE')
  call plot1d(50,51,1,rhomt,rhoir)
  close(50)
  close(51)
end if
! charge density in 2D
if (tdrho2d) then
  open(50,file='RHO2D'//trim(fext),form='FORMATTED',action='WRITE')
  call plot2d(.false.,50,1,rhomt,rhoir)
  close(50)
end if
! charge density in 3D
if (tdrho3d) then
  open(50,file='RHO3D'//trim(fext),form='FORMATTED',action='WRITE')
  call plot3d(50,1,rhomt,rhoir)
  close(50)
end if
! magnetisation in 1D, 2D or 3D
if ((tdmag1d.or.tdmag2d.or.tdmag3d).and.(spinpol)) then
  allocate(rvfmt(npmtmax,natmtot,3),rvfir(ngtot,3))
  if (ncmag) then
! non-collinear
    rvfmt(:,:,:)=magmt(:,:,:)
    rvfir(:,:)=magir(:,:)
  else
! collinear
    rvfmt(:,:,1:2)=0.d0
    rvfir(:,1:2)=0.d0
    rvfmt(:,:,3)=magmt(:,:,1)
    rvfir(:,3)=magir(:,1)
  end if
  if (tdmag1d) then
    open(50,file='MAG1D'//trim(fext),form='FORMATTED',action='WRITE')
    open(51,file='MAGLINES.OUT',form='FORMATTED',action='WRITE')
    call plot1d(50,51,3,rvfmt,rvfir)
    close(50)
    close(51)
  end if
  if (tdmag2d) then
    open(50,file='MAG2D'//trim(fext),form='FORMATTED',action='WRITE')
    call plot2d(.true.,50,3,rvfmt,rvfir)
    close(50)
  end if
  if (tdmag3d) then
    open(50,file='MAG3D'//trim(fext),form='FORMATTED',action='WRITE')
    call plot3d(50,3,rvfmt,rvfir)
    close(50)
  end if
  deallocate(rvfmt,rvfir)
end if
! gauge-invariant current density in 1D
if (tdjr1d) then
  open(50,file='JR1D'//trim(fext),form='FORMATTED',action='WRITE')
  open(51,file='JRLINES.OUT',form='FORMATTED',action='WRITE')
  call plot1d(50,51,3,jrmt,jrir)
  close(50)
  close(51)
end if
! gauge-invariant current density in 2D
if (tdjr2d) then
  open(50,file='JR2D'//trim(fext),form='FORMATTED',action='WRITE')
  call plot2d(.true.,50,3,jrmt,jrir)
  close(50)
end if
! gauge-invariant current density in 3D
if (tdjr3d) then
  open(50,file='JR3D'//trim(fext),form='FORMATTED',action='WRITE')
  call plot3d(50,3,jrmt,jrir)
  close(50)
end if
! calculate and write tensor moments
if (dftu /= 0) then
  if (tmwrite) call writetdtm3
end if
end subroutine

