
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepmain
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,idm,is,ias
integer nrc,nrci,np,npc
integer n,nthd,it
real(8) resp,t1
! automatic arrays
real(8) rfmt2(npcmtmax)
! allocatable arrays
real(8), allocatable :: rfmt1(:,:),rfir(:)
real(8), allocatable :: rvfmt(:,:,:),rvfir(:,:)
complex(8), allocatable :: vclcv(:,:,:,:),vclvv(:,:,:)
! external functions
real(8), external :: rfint,rfinpc
if (iscl <= 0) return
! calculate Coulomb matrix elements
allocate(vclcv(ncrmax,natmtot,nstsv,nkpt),vclvv(nstsv,nstsv,nkpt))
call oepvcl(vclcv,vclvv)
! allocate local arrays
allocate(rfmt1(npmtmax,natmtot),rfir(ngtot))
if (spinpol) then
  allocate(rvfmt(npmtmax,natmtot,ndmag),rvfir(ngtot,ndmag))
end if
! zero initial exchange potential and magnetic field
vxmt(:,:)=0.d0
vxir(:)=0.d0
if (spinpol) then
  bxmt(:,:,:)=0.d0
  bxir(:,:)=0.d0
end if
!------------------------------!
!     start iteration loop     !
!------------------------------!
do it=1,maxitoep
  if (mp_mpi.and.(mod(it,10) == 0)) then
    write(*,'("Info(oepmain): done ",I4," iterations of ",I4)') it,maxitoep
  end if
! zero the residuals
  dvxmt(:,:)=0.d0
  dvxir(:)=0.d0
  if (spinpol) then
    dbxmt(:,:,:)=0.d0
    dbxir(:,:)=0.d0
  end if
! calculate the k-dependent residuals
  call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
  do ik=1,nkpt
! distribute among MPI processes
    if (mod(ik-1,np_mpi) /= lp_mpi) cycle
    call oepresk(ik,vclcv,vclvv)
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
! add residuals from each process and redistribute
  if (np_mpi > 1) then
    n=npcmtmax*natmtot
    call mpi_allreduce(mpi_in_place,dvxmt,n,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
    call mpi_allreduce(mpi_in_place,dvxir,ngtot,mpi_double_precision, &
     mpi_sum,mpicom,ierror)
    if (spinpol) then
      n=n*ndmag
      call mpi_allreduce(mpi_in_place,dbxmt,n,mpi_double_precision,mpi_sum, &
       mpicom,ierror)
      n=ngtot*ndmag
      call mpi_allreduce(mpi_in_place,dbxir,n,mpi_double_precision,mpi_sum, &
       mpicom,ierror)
    end if
  end if
! convert muffin-tin residuals to spherical harmonics
  call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is,nrc,nrci,idm) &
!$OMP NUM_THREADS(nthd)
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    call rfsht(nrc,nrci,dvxmt(:,ias),rfmt1(:,ias))
    do idm=1,ndmag
      call rfsht(nrc,nrci,dbxmt(:,ias,idm),rvfmt(:,ias,idm))
    end do
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
! symmetrise the residuals
  call symrf(nrcmt,nrcmti,npcmt,ngridg,ngtot,ngvec,igfft,npmtmax,rfmt1,dvxir)
  if (spinpol) then
    call symrvf(.true.,ncmag,nrcmt,nrcmti,npcmt,ngridg,ngtot,ngvec,igfft, &
     npmtmax,rvfmt,ngtot,dbxir)
  end if
! magnitude of residuals
  resoep=rfinpc(npmtmax,rfmt1,dvxir,rfmt1,dvxir)
  do idm=1,ndmag
    t1=rfinpc(npmtmax,rvfmt(:,:,idm),dbxir(:,idm),rvfmt(:,:,idm),dbxir(:,idm))
    resoep=resoep+t1
  end do
  resoep=sqrt(resoep)/omega
! update the step size
  if (it >= 2) then
! check residual against previous value
    if (resoep < resp) then
      tauoep=tauoep+tau0oep
    else
      tauoep=tau0oep
    end if
  end if
! store previous residual
  resp=resoep
! update exchange potential and magnetic field
  call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt2,is,nrc,nrci,npc,idm) &
!$OMP NUM_THREADS(nthd)
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
! convert residual to spherical coordinates
    call rbsht(nrc,nrci,rfmt1(:,ias),rfmt2)
! subtract from exchange potential
    vxmt(1:npc,ias)=vxmt(1:npc,ias)-tauoep*rfmt2(1:npc)
! repeat for exchange magnetic field
    do idm=1,ndmag
      call rbsht(nrc,nrci,rvfmt(:,ias,idm),rfmt2)
      bxmt(1:npc,ias,idm)=bxmt(1:npc,ias,idm)-tauoep*rfmt2(1:npc)
    end do
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
  vxir(:)=vxir(:)-tauoep*dvxir(:)
  do idm=1,ndmag
    bxir(:,idm)=bxir(:,idm)-tauoep*dbxir(:,idm)
  end do
! end iteration loop
end do
! convert the exchange potential and field to spherical harmonics
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is,nrc,nrci,idm) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  call rfsht(nrc,nrci,vxmt(:,ias),rfmt1(:,ias))
  do idm=1,ndmag
    call rfsht(nrc,nrci,bxmt(:,ias,idm),rvfmt(:,ias,idm))
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! convert potential and field from a coarse to a fine radial mesh
call rfmtctof(rfmt1)
do idm=1,ndmag
  call rfmtctof(rvfmt(:,:,idm))
end do
! add to existing (density derived) correlation potential and field
do ias=1,natmtot
  is=idxis(ias)
  np=npmt(is)
  vxcmt(1:np,ias)=vxcmt(1:np,ias)+rfmt1(1:np,ias)
  do idm=1,ndmag
    bxcmt(1:np,ias,idm)=bxcmt(1:np,ias,idm)+rvfmt(1:np,ias,idm)
  end do
end do
vxcir(:)=vxcir(:)+vxir(:)
do idm=1,ndmag
  bxcir(:,idm)=bxcir(:,idm)+bxir(:,idm)
end do
! symmetrise the exchange potential and field
call symrf(nrmt,nrmti,npmt,ngridg,ngtot,ngvec,igfft,npmtmax,vxcmt,vxcir)
if (spinpol) then
  call symrvf(.true.,ncmag,nrmt,nrmti,npmt,ngridg,ngtot,ngvec,igfft,npmtmax, &
   bxcmt,ngtot,bxcir)
end if
deallocate(rfmt1,rfir,vclcv,vclvv)
if (spinpol) deallocate(rvfmt,rvfir)
! set the constant part of the exchange potential equal to zero
call rfint0(0.d0,vxcmt,vxcir)
end subroutine

