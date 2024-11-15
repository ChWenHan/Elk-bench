
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetdjtk
use modmain
use modtddft
use modmpi
use modomp
implicit none
! local variables
integer ik,ist,l,lp,nthd
real(8) ca,t1,t2
complex(8) z1
character(32) fext
! allocatable arrays
real(8), allocatable :: jtk(:,:)
complex(8), allocatable :: evecsv(:,:),evecsvt(:,:)
complex(8), allocatable :: pmat(:,:,:),a(:,:),b(:,:)
! external functions
complex(8), external :: zdotc
! coupling constant of the external A-field (-1/c)
ca=-1.d0/solsc
allocate(jtk(3,nkpt))
jtk(:,:)=0.d0
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecsv,evecsvt,pmat) &
!$OMP PRIVATE(a,b,l,t1,t2,ist,z1) &
!$OMP NUM_THREADS(nthd)
allocate(evecsv(nstsv,nstsv),evecsvt(nstsv,nstsv))
allocate(pmat(nstsv,nstsv,3),a(nstsv,nstsv),b(nstsv,nstsv))
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
! get the momentum matrix elements from file
  call getpmat(vkl(:,ik),pmat)
! read in ground-state eigenvectors
  call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! read in time-dependent Kohn-Sham eigenvectors (first-variational basis)
  call getevecsv(filext,ik,vkl(:,ik),evecsvt)
  do l=1,3
! form the momentum operator matrix elements in the first-variational basis
    call zgemm('N','C',nstsv,nstsv,nstsv,zone,pmat(:,:,l),nstsv,evecsv,nstsv, &
     zzero,a,nstsv)
    call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,a,nstsv,zzero,b, &
     nstsv)
    call zgemm('N','N',nstsv,nstsv,nstsv,zone,b,nstsv,evecsvt,nstsv,zzero,a, &
     nstsv)
! add to the total current for this k-point (including diamagnetic contribution)
    t1=ca*afieldt(l,itimes)
    do ist=1,nstsv
      t2=occsv(ist,ik)
      if (abs(t2) > epsocc) then
        z1=zdotc(nstsv,evecsvt(:,ist),1,a(:,ist),1)
        jtk(l,ik)=jtk(l,ik)+t2*(dble(z1)+t1)
      end if
    end do
  end do
end do
!$OMP END DO
deallocate(evecsv,evecsvt,pmat,a,b)
!$OMP END PARALLEL
call freethd(nthd)
! broadcast current array to every MPI process
if (np_mpi > 1) then
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(jtk(:,ik),3,mpi_double_precision,lp,mpicom,ierror)
  end do
end if
! write k-point dependent total current to file
if (mp_mpi) then
! file extension
  write(fext,'("_TS",I8.8,".OUT")') itimes
  open(50,file='JTOTK'//trim(fext),form='FORMATTED',action='WRITE')
  do ik=1,nkpt
    write(50,'(I6,6G18.10)') ik,vkl(:,ik),jtk(:,ik)
  end do
  close(50)
end if
deallocate(jtk)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

