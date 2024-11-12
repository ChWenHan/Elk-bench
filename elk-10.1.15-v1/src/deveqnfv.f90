
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine deveqnfv(ngp,ngpq,igpig,igpqig,vgpc,vgpqc,evalfv,apwalm,apwalmq, &
 dapwalm,dapwalmq,evecfv,devalfvp,devecfv)
use modmain
use modphonon
use modomp
implicit none
! arguments
integer, intent(in) :: ngp,ngpq
integer, intent(in) :: igpig(ngkmax),igpqig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax),vgpqc(3,ngkmax)
real(8), intent(in) :: evalfv(nstfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalmq(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
real(8), intent(out) :: devalfvp(nstfv)
complex(8), intent(out) :: devecfv(nmatmax,nstfv)
! local variables
integer nm,nmq,is,ias,jst,i
integer lwork,info,nthd
real(8) t1
complex(8) z1
! allocatable arrays
real(8), allocatable :: w(:),rwork(:)
complex(8), allocatable :: h(:,:),o(:,:),dh(:,:),od(:,:)
complex(8), allocatable :: x(:),y(:),work(:)
! external functions
real(8), external :: ddot
! matrix sizes for k and k+q
nm=ngp+nlotot
nmq=ngpq+nlotot
allocate(h(nmq,nmq),o(nmq,nmq))
! compute the Hamiltonian and overlap matrices at p+q
call holdthd(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
call hmlfv(nmq,ngpq,igpqig,vgpqc,apwalmq,h)
!$OMP SECTION
call olpfv(nmq,ngpq,igpqig,apwalmq,o)
!$OMP END PARALLEL SECTIONS
call freethd(nthd)
! solve the generalised eigenvalue problem (H - eⱼ O)|vⱼ> = 0
! (note: these are also the eigenvalues/vectors of O⁻¹ H )
lwork=2*nmq
allocate(w(nmq),rwork(3*nmq),work(lwork))
call zhegv(1,'V','U',nmq,h,nmq,o,nmq,w,work,lwork,rwork,info)
if (info /= 0) then
  write(*,*)
  write(*,'("Error(deveqnfv): diagonalisation failed")')
  write(*,'(" ZHEGV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
deallocate(rwork,o,work)
! compute the Hamiltonian and overlap matrix derivatives
allocate(dh(nmq,nm),od(nmq,nm))
call holdthd(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP PRIVATE(ias,is) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
call dhmlistl(ngp,ngpq,igpig,igpqig,vgpc,vgpqc,nmq,dh)
dh(ngpq+1:nmq,1:ngp)=0.d0
dh(1:nmq,ngp+1:nm)=0.d0
do ias=1,natmtot
  is=idxis(ias)
  call dhmlaa(is,ias,ngp,ngpq,apwalm(:,:,:,ias),apwalmq(:,:,:,ias),dapwalm, &
   dapwalmq,nmq,dh)
  call dhmlalo(is,ias,ngp,ngpq,apwalm(:,:,:,ias),apwalmq(:,:,:,ias),dapwalm, &
   dapwalmq,nmq,dh)
  call dhmllolo(is,ias,ngp,ngpq,nmq,dh)
end do
!$OMP SECTION
call dolpistl(ngp,ngpq,igpig,igpqig,nmq,od)
od(ngpq+1:nmq,1:ngp)=0.d0
od(1:nmq,ngp+1:nm)=0.d0
do ias=1,natmtot
  is=idxis(ias)
  call dolpaa(is,ias,ngp,ngpq,apwalm(:,:,:,ias),apwalmq(:,:,:,ias),dapwalm, &
   dapwalmq,nmq,od)
  call dolpalo(is,ias,ngp,ngpq,dapwalm,dapwalmq,nmq,od)
end do
!$OMP END PARALLEL SECTIONS
call freethd(nthd)
allocate(x(nmq),y(nmq))
! loop over states
do jst=1,nstfv
! compute |dvⱼ> = V (eⱼ - D)⁻¹ V† (dH - eⱼ dO)|vⱼ>
  z1=-evalfv(jst)
  call zgemv('N',nmq,nm,z1,od,nmq,evecfv(:,jst),1,zzero,x,1)
  call zgemv('N',nmq,nm,zone,dh,nmq,evecfv(:,jst),1,zone,x,1)
! compute the first-order change in eigenvalue
  if (tphq0) then
    devalfvp(jst)=ddot(2*nmq,evecfv(:,jst),1,x,1)
  else
    devalfvp(jst)=0.d0
  end if
  call zgemv('C',nmq,nmq,zone,h,nmq,x,1,zzero,y,1)
  do i=1,nmq
    t1=evalfv(jst)-w(i)
    if (abs(t1) > epsdev) then
      y(i)=y(i)/t1
    else
      y(i)=0.d0
    end if
  end do
  call zgemv('N',nmq,nmq,zone,h,nmq,y,1,zzero,devecfv(:,jst),1)
end do
deallocate(w,h,dh,od,x,y)
end subroutine

