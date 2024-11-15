
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dforce(dyn)
use modmain
use modphonon
use modmpi
use modomp
implicit none
! arguments
complex(8), intent(out) :: dyn(3,natmtot)
! local variables
integer ik,is,ias,jas,ip
integer nr,nri,nthd
complex(8) z1
! automatic arrays
complex(8) dynibs(3,natmtot)
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: grhomt(:,:,:),grhoir(:,:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
complex(8), allocatable :: gvclmt(:,:,:),gvclir(:,:)
complex(8), allocatable :: zfmt(:),gvcln(:,:),gzfmt(:,:)
! external functions
complex(8), external :: zfmtinp
allocate(zrhomt(npmtmax,natmtot),zrhoir(ngtot))
allocate(grhomt(npmtmax,natmtot,3),grhoir(ngtot,3))
allocate(zvclmt(npmtmax,natmtot+1),zvclir(ngtot))
allocate(gvclmt(npmtmax,natmtot,3),gvclir(ngtot,3))
allocate(zfmt(npmtmax),gvcln(npmtmax,3),gzfmt(npmtmax,3))
! make complex copy of the density
do ias=1,natmtot
  is=idxis(ias)
  call rtozfmt(nrmt(is),nrmti(is),rhomt(:,ias),zrhomt(:,ias))
end do
zrhoir(1:ngtot)=rhoir(1:ngtot)
! compute the gradient of the density
call gradzf(zrhomt,zrhoir,grhomt,grhoir)
!--------------------------------------------------------------!
!     Hellmann-Feynman force derivative for displaced atom     !
!--------------------------------------------------------------!
! calculate the gradient of the nuclear potential
call gradzvcln(isph,gvcln)
do ip=1,3
! compute the q-dependent nuclear Coulomb potential derivative
  zvclmt(1:npmtmax,1:natmtot)=0.d0
  zvclmt(1:npmt(isph),iasph)=gvcln(1:npmt(isph),ip)
! zero the interstitial density
  zrhoir(1:ngtot)=0.d0
  if (ip == ipph) then; jas=iasph; else; jas=0; endif
  call zpotcoul(jas,nrmt,nrmti,npmt,nrmtmax,rlmt,ngridg,igfft,ngvec,gqc,gclgq, &
   ngvec,jlgqrmt,ylmgq,sfacgq,zrhoir,npmtmax,zvclmt,zvclir)
! multiply with density derivative and integrate
  z1=sum(cfunir(1:ngtot)*conjg(zvclir(1:ngtot))*drhoir(1:ngtot))
  z1=z1*omega/dble(ngtot)
  do ias=1,natmtot
    is=idxis(ias)
    z1=z1+zfmtinp(nrmt(is),nrmti(is),wrmt(:,is),zvclmt(:,ias),drhomt(:,ias))
  end do
  dyn(ip,iasph)=-z1
! nuclear-nuclear term
  if (ip == ipph) then; jas=natmtot+1; else; jas=iasph; end if
  call gradzfmt(nrmt(isph),nrmti(isph),rlmt(:,-1,isph),wcrmt(:,:,isph), &
   zvclmt(:,jas),npmtmax,gzfmt)
  z1=spzn(isph)*gzfmt(1,ipph)*y00
  dyn(ip,iasph)=dyn(ip,iasph)-z1
end do
! compute the lattice-periodic nuclear Coulomb potential derivative
zvclmt(1:npmtmax,1:natmtot)=0.d0
zvclmt(1:npmt(isph),iasph)=gvcln(1:npmt(isph),ipph)
call zpotcoul(iasph,nrmt,nrmti,npmt,nrmtmax,rlmt,ngridg,igfft,ngvec,gc,gclg, &
 ngvec,jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)
do ip=1,3
! multiply with density gradient and integrate
  z1=sum(cfunir(1:ngtot)*zvclir(1:ngtot)*grhoir(1:ngtot,ip))
  z1=z1*omega/dble(ngtot)
  do ias=1,natmtot
    is=idxis(ias)
    z1=z1+zfmtinp(nrmt(is),nrmti(is),wrmt(:,is),grhomt(:,ias,ip),zvclmt(:,ias))
  end do
  dyn(ip,iasph)=dyn(ip,iasph)-z1
end do
! nuclear-nuclear term
do ip=1,3
  do ias=1,natmtot
    is=idxis(ias)
    if (ias == iasph) then; jas=natmtot+1; else; jas=ias; end if
    call gradzfmt(nrmt(is),nrmti(is),rlmt(:,-1,is),wcrmt(:,:,is),zvclmt(:,jas),&
     npmtmax,gzfmt)
    z1=spzn(is)*gzfmt(1,ip)*y00
    dyn(ip,iasph)=dyn(ip,iasph)+z1
  end do
end do
!-------------------------------------------------------------------!
!     Hellmann-Feynman force derivative for non-displaced atoms     !
!-------------------------------------------------------------------!
do ias=1,natmtot
  if (ias == iasph) cycle
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
! compute the gradient of the Coulomb potential derivative at the nucleus
  call gradzfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),dvclmt(:,ias),npmtmax,gzfmt)
  do ip=1,3
    dyn(ip,ias)=spzn(is)*gzfmt(1,ip)*y00
  end do
end do
!--------------------------------------------!
!     IBS correction to force derivative     !
!--------------------------------------------!
dynibs(:,:)=0.d0
! k-point dependent part
call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP REDUCTION(+:dynibs) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
  call dforcek(ik,dynibs)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! add IBS force derivatives from each process and redistribute
if (np_mpi > 1) then
  call mpi_allreduce(mpi_in_place,dynibs,3*natmtot,mpi_double_complex,mpi_sum, &
   mpicom,ierror)
end if
! k-point independent part
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  do ip=1,3
    z1=zfmtinp(nr,nri,wrmt(:,is),grhomt(:,ias,ip),dvsmt(:,ias))
    dynibs(ip,ias)=dynibs(ip,ias)-z1
  end do
! compute the gradient of the density derivative
  if (ias == iasph) then; jas=natmtot+1; else; jas=ias; end if
  call gradzfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),drhomt(:,jas),npmtmax,gzfmt)
! convert Kohn-Sham potential to complex spherical harmonics
  call rtozfmt(nr,nri,vsmt(:,ias),zfmt)
  do ip=1,3
    z1=zfmtinp(nr,nri,wrmt(:,is),zfmt,gzfmt(:,ip))
    dynibs(ip,ias)=dynibs(ip,ias)-z1
  end do
end do
! add the IBS force derivatives to the total
dyn(:,:)=dyn(:,:)+dynibs(:,:)
deallocate(zrhomt,zrhoir,grhomt,grhoir)
deallocate(zvclmt,zvclir,gvclmt,gvclir,zfmt,gzfmt)
end subroutine

