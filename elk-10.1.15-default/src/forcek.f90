
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: forcek
! !INTERFACE:
subroutine forcek(ik,forceibs)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   ik       : reduced k-point number (in,integer)
!   forceibs : IBS force (inout,real(3,natmtot))
! !DESCRIPTION:
!   Computes the {\bf k}-dependent contribution to the incomplete basis set
!   (IBS) force. See the calling routine {\tt force} for a full description.
!
! !REVISION HISTORY:
!   Created June 2006 (JKD)
!   Updated for spin-spiral case, May 2007 (Francesco Cricchio and JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(inout) :: forceibs(3,natmtot)
! local variables
integer ispn0,ispn1,ispn,jspn
integer n,nm,is,ias,ist,jst
integer j1,j2,j3,ig,i,j,l,nthd
real(8) v1,v2,v3,sm,t1
complex(8) z1,z2
! automatic arrays
real(8) evalfv(nstfv,nspnfv)
complex(8) vh(nmatmax),vo(nmatmax)
complex(8) ffv(nstfv,nstfv),y(nstfv)
! allocatable arrays
integer, allocatable :: ijg(:,:)
real(8), allocatable :: dp(:,:)
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: h(:,:),o(:,:),dlh(:,:),dlo(:,:)
! external functions
complex(8), external :: zdotc
! allocate local arrays
allocate(ijg(nmatmax,nmatmax),dp(nmatmax,nmatmax))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(h(nmatmax,nmatmax),o(nmatmax,nmatmax))
allocate(dlh(nmatmax,nmatmax),dlo(nmatmax,nmatmax))
! get the eigenvalues/vectors from file
call getevalfv(filext,ik,vkl(:,ik),evalfv)
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
if (tevecsv) then
  allocate(evecsv(nstsv,nstsv))
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
end if
! loop over first-variational spin components
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
  n=ngk(jspn,ik)
  nm=nmat(jspn,ik)
  do j=1,n
    ig=igkig(j,jspn,ik)
    j1=ivg(1,ig); j2=ivg(2,ig); j3=ivg(3,ig)
    v1=0.5d0*vgkc(1,j,jspn,ik)
    v2=0.5d0*vgkc(2,j,jspn,ik)
    v3=0.5d0*vgkc(3,j,jspn,ik)
    do i=1,j
      ig=igkig(i,jspn,ik)
      ijg(i,j)=ivgig(ivg(1,ig)-j1,ivg(2,ig)-j2,ivg(3,ig)-j3)
      dp(i,j)=vgkc(1,i,jspn,ik)*v1+vgkc(2,i,jspn,ik)*v2+vgkc(3,i,jspn,ik)*v3
    end do
  end do
! find the matching coefficients
  call match(n,vgkc(:,:,jspn,ik),gkc(:,jspn,ik),sfacgk(:,:,jspn,ik),apwalm)
! zero the local-orbital-local-orbital contribution
  do j=n+1,nm
    dlh(n+1:j,j)=0.d0
    dlo(n+1:j,j)=0.d0
  end do
! loop over species and atoms
  do ias=1,natmtot
    is=idxis(ias)
! Hamiltonian and overlap matrices
    do j=1,nm
      h(1:j,j)=0.d0
    end do
    call hmlaa(.false.,is,ias,n,apwalm(:,:,:,ias),nmatmax,h)
    call hmlalo(is,ias,n,apwalm(:,:,:,ias),nmatmax,h)
    do j=1,nm
      o(1:j,j)=0.d0
    end do
    call olpaa(.false.,is,n,apwalm(:,:,:,ias),nmatmax,o)
    call olpalo(is,ias,n,apwalm(:,:,:,ias),nmatmax,o)
! loop over Cartesian directions
    do l=1,3
! APW-APW contribution
      do j=1,n
        do i=1,j
          ig=ijg(i,j)
          t1=vgc(l,ig)
          z1=-ffacg(ig,is)*conjg(sfacg(ig,ias))
          z2=t1*(dp(i,j)*z1+h(i,j))
          dlh(i,j)=cmplx(-aimag(z2),dble(z2),8)
          z2=t1*(z1+o(i,j))
          dlo(i,j)=cmplx(-aimag(z2),dble(z2),8)
        end do
      end do
! APW-local-orbital contribution
      do j=n+1,nm
        do i=1,n
          t1=vgkc(l,i,jspn,ik)
          z1=t1*h(i,j)
          dlh(i,j)=cmplx(-aimag(z1),dble(z1),8)
          z1=t1*o(i,j)
          dlo(i,j)=cmplx(-aimag(z1),dble(z1),8)
        end do
      end do
! compute the force matrix elements in the first-variational basis
      call holdthd(nstfv,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(vh,vo,t1,ist,z1,z2) &
!$OMP NUM_THREADS(nthd)
      do jst=1,nstfv
        call zhemv('U',nm,zone,dlh,nmatmax,evecfv(:,jst,jspn),1,zzero,vh,1)
        call zhemv('U',nm,zone,dlo,nmatmax,evecfv(:,jst,jspn),1,zzero,vo,1)
        t1=evalfv(jst,jspn)
        do ist=1,nstfv
          z1=zdotc(nm,evecfv(:,ist,jspn),1,vh,1)
          z2=zdotc(nm,evecfv(:,ist,jspn),1,vo,1)
          ffv(ist,jst)=z1-t1*z2
        end do
      end do
!$OMP END PARALLEL DO
      call freethd(nthd)
! compute the force using the second-variational coefficients if required
      sm=0.d0
      if (tevecsv) then
! spin-polarised case
        do j=1,nstsv
          do ispn=ispn0,ispn1
            i=(ispn-1)*nstfv+1
            call zgemv('N',nstfv,nstfv,zone,ffv,nstfv,evecsv(i,j),1,zzero,y,1)
            z1=zdotc(nstfv,evecsv(i,j),1,y,1)
            sm=sm+occsv(j,ik)*dble(z1)
          end do
        end do
      else
! spin-unpolarised case
        do j=1,nstsv
          sm=sm+occsv(j,ik)*dble(ffv(j,j))
        end do
      end if
      forceibs(l,ias)=forceibs(l,ias)+wkpt(ik)*sm
! end loop over Cartesian components
    end do
! end loop over atoms and species
  end do
! end loop over first-variational spins
end do
deallocate(ijg,dp,apwalm,evecfv)
deallocate(h,o,dlh,dlo)
if (tevecsv) deallocate(evecsv)
end subroutine
!EOC

