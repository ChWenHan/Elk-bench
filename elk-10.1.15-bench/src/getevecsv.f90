
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevecsv(fext,ikp,vpl,evecsv)
use modmain
use modramdisk
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ikp
real(8), intent(in) :: vpl(3)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
logical tgs
integer isym,lspn
integer ik,ist,jst
integer recl,nstsv_
real(8) vkl_(3),det,v(3),th,t1
complex(8) su2(2,2),z1,z2
character(256) fname
if (ikp > 0) then
  ik=ikp
else
! find the equivalent k-point number and symmetry which rotates vkl to vpl
  call findkpt(vpl,isym,ik)
end if
! construct the filename
fname=trim(scrpath)//'EVECSV'//trim(fext)
!$OMP CRITICAL(u206)
! read from RAM disk if required
if (ramdisk) then
  call getrd(fname,ik,tgs,v1=vkl_,n1=nstsv_,nzv=nstsv*nstsv,zva=evecsv)
  if (tgs) goto 10
end if
! find the record length
inquire(iolength=recl) vkl_,nstsv_,evecsv
open(206,file=fname,form='UNFORMATTED',access='DIRECT',recl=recl)
read(206,rec=ik) vkl_,nstsv_,evecsv
close(206)
10 continue
!$OMP END CRITICAL(u206)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1 > epslat) then
  write(*,*)
  write(*,'("Error(getevecsv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVECSV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv /= nstsv_) then
  write(*,*)
  write(*,'("Error(getevecsv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" EVECSV.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
! if eigenvectors are spin-unpolarised return
if (.not.spinpol) return
! if p = k then return
if (ikp > 0) return
! index to global spin rotation in lattice point group
lspn=lspnsymc(isym)
! if symmetry element is the identity return
if (lspn == 1) return
! find the SU(2) representation of the spin rotation matrix
call rotaxang(epslat,symlatc(:,:,lspn),det,v,th)
call axangsu2(v,th,su2)
! apply SU(2) matrix to second-variational states (active transformation)
do jst=1,nstsv
  do ist=1,nstfv
    z1=evecsv(ist,jst)
    z2=evecsv(ist+nstfv,jst)
    evecsv(ist,jst)=su2(1,1)*z1+su2(1,2)*z2
    evecsv(ist+nstfv,jst)=su2(2,1)*z1+su2(2,2)*z2
  end do
end do
end subroutine

