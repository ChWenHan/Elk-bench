
! Copyright (C) 2024 Eddie Harris-Lee, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine energyafsp
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,ist,i,l,nthd
real(8) ca,sm,v(3),wo,t1
complex(8) z11,z12,z21,z22
! automatic arrays
complex(8) y(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgk(:,:,:),pmat(:,:,:)
! external functions
complex(8), external :: zdotc
! coupling constant of the external spin-polarised A-field (-1/c)
ca=-1.d0/solsc
sm=0.d0
! begin parallel loop over k-points
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv) &
!$OMP PRIVATE(wfmt,wfgk,pmat,y) &
!$OMP PRIVATE(ist,wo,l,v,z11,z12,z21,z22,i,t1) &
!$OMP REDUCTION(+:sm) &
!$OMP NUM_THREADS(nthd)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfgk(ngkmax,nspinor,nstsv))
allocate(pmat(nstsv,nstsv,3))
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
! get the eigenvectors from file
  call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
! find the matching coefficients
  call match(ngk(1,ik),vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
! calculate the wavefunctions for all states
  call genwfsv(.true.,.true.,nstsv,[0],ngridg,igfft,ngk(:,ik),igkig(:,:,ik), &
   apwalm,evecfv,evecsv,wfmt,ngkmax,wfgk)
! calculate the momentum matrix elements in the second-variational basis
  call genpmatk(ngk(:,ik),igkig(:,:,ik),vgkc(:,:,:,ik),wfmt,wfgk,pmat)
  do ist=1,nstsv
    wo=occsv(ist,ik)
    if (abs(wo) < epsocc) cycle
    wo=wo*wkpt(ik)
    do l=1,3
      v(1:3)=ca*afspc(l,1:3)
      z11=v(3)
      z12=cmplx(v(1),-v(2),8)
      z21=cmplx(v(1),v(2),8)
      z22=-v(3)
! convert momentum matrix elements to first-variational basis
      call zgemv('N',nstsv,nstsv,zone,evecsv,nstsv,pmat(:,ist,l),1,zzero,y,1)
      i=nstfv+1
      t1=dble(z11*zdotc(nstfv,evecsv(1,ist),1,y,1)) &
        +dble(z12*zdotc(nstfv,evecsv(1,ist),1,y(i),1)) &
        +dble(z21*zdotc(nstfv,evecsv(i,ist),1,y,1)) &
        +dble(z22*zdotc(nstfv,evecsv(i,ist),1,y(i),1))
! add to the expectation value
      sm=sm+wo*t1
    end do
  end do
end do
!$OMP END DO
deallocate(apwalm,evecfv,evecsv)
deallocate(wfmt,wfgk,pmat)
!$OMP END PARALLEL
call freethd(nthd)
! add sums from each process and redistribute
if (np_mpi > 1) then
  call mpi_allreduce(mpi_in_place,sm,1,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
end if
! subtract from kinetic energy
engykn=engykn-0.5d0*sm
end subroutine

