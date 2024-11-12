
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gencore
! !INTERFACE:
subroutine gencore
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Computes the core radial wavefunctions, eigenvalues and densities. The
!   radial Dirac equation is solved in the spherical part of the Kohn-Sham
!   potential to which the atomic potential has been appended for
!   $r>R_{\rm MT}$. In the case of spin-polarised calculations, and when
!   {\tt spincore} is set to {\tt .true.}, the Dirac equation is solved in the
!   spin-up and -down potentials created from the Kohn-Sham scalar potential and
!   magnetic field magnitude, with the occupancy divided equally between up and
!   down. The up and down densities determined in this way are added to both the
!   scalar density and the magnetisation in the routine {\tt rhocore}. Note
!   that this procedure is a simple, but inexact, approach to solving the radial
!   Dirac equation in a magnetic field.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Added polarised cores, November 2009 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ist,ispn,idm
integer is,ia,ja,ias,jas
integer nr,nri,nrs,nthd
real(8) v(ndmag),t1
! automatic arrays
logical done(natmtot)
real(8) vr(nrspmax),br(nrmtmax),fr(nrmtmax)
real(8) eval(nstspmax)
! external functions
real(8), external :: rfmtint
done(:)=.false.
! loop over species and atoms
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(vr,br,fr,eval) &
!$OMP PRIVATE(is,ia,ja,jas) &
!$OMP PRIVATE(nr,nri,nrs,idm) &
!$OMP PRIVATE(v,t1,ispn,ist) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
! flag all equivalent atoms as done
!$OMP CRITICAL
  do ja=1,natoms(is)
    jas=idxas(ja,is)
    if ((ia /= ja).and.(eqatoms(ia,ja,is))) done(jas)=.true.
  end do
!$OMP END CRITICAL
  if (done(ias)) cycle
  nr=nrmt(is)
  nri=nrmti(is)
  nrs=nrsp(is)
! Kohn-Sham magnetic field for spin-polarised core
  if (spincore) then
! compute the averaged direction of the exchange-correlation magnetic field
    do idm=1,ndmag
      v(idm)=rfmtint(nr,nri,wrmt(:,is),bxcmt(:,ias,idm))
    end do
! normalise
    if (ncmag) then
      t1=sqrt(v(1)**2+v(2)**2+v(3)**2)
    else
      t1=abs(v(1))
    end if
    if (t1 > 1.d-10) v(1:ndmag)=v(1:ndmag)/t1
! determine the component of the field along the averaged direction
    br(1:nr)=0.d0
    do idm=1,ndmag
! extract the spherical (l = m = 0) component of B_xc
      call rfmtlm(1,nr,nri,bxcmt(:,ias,idm),fr)
      t1=v(idm)*y00
      br(1:nr)=br(1:nr)+t1*fr(1:nr)
    end do
  end if
! loop over spin channels
  do ispn=1,nspncr
! use the spherical part of the crystal Kohn-Sham potential
    call rfmtlm(1,nr,nri,vsmt(:,ias),vr)
    vr(1:nr)=vr(1:nr)*y00
! spin-up and -down potentials for polarised core
    if (spincore) then
      if (ispn == 1) then
        vr(1:nr)=vr(1:nr)-br(1:nr)
      else
        vr(1:nr)=vr(1:nr)+br(1:nr)
      end if
    end if
! append the Kohn-Sham potential from the atomic calculation for r > R_MT
    t1=vr(nr)-vrsp(nr,is)
    vr(nr+1:nrs)=vrsp(nr+1:nrs,is)+t1
    rhocr(1:nr,ias,ispn)=0.d0
    do ist=1,nstsp(is)
      if (spcore(ist,is)) then
! solve the Dirac equation
        eval(ist)=evalcr(ist,ias)
        call rdirac(solsc,nsp(ist,is),lsp(ist,is),ksp(ist,is),nrs,rsp(:,is),vr,&
         eval(ist),rwfcr(:,1,ist,ias),rwfcr(:,2,ist,ias))
        if (spincore) then
! use the spin-averaged eigenvalue for the polarised core
          if (ispn == 1) then
            evalcr(ist,ias)=eval(ist)
          else
            evalcr(ist,ias)=0.5d0*(evalcr(ist,ias)+eval(ist))
          end if
          t1=0.5d0*occcr(ist,ias)
        else
          evalcr(ist,ias)=eval(ist)
          t1=occcr(ist,ias)
        end if
! add to the core density
        rhocr(1:nr,ias,ispn)=rhocr(1:nr,ias,ispn) &
         +t1*(rwfcr(1:nr,1,ist,ias)**2+rwfcr(1:nr,2,ist,ias)**2)
      end if
    end do
    rhocr(1:nr,ias,ispn)=rhocr(1:nr,ias,ispn)*rlmt(1:nr,-2,is)*y00
! end loop over spin channels
  end do
! flag atom as done and copy to equivalent atoms
  done(ias)=.true.
  do ja=1,natoms(is)
    jas=idxas(ja,is)
    if (eqatoms(ia,ja,is)) then
      do ist=1,nstsp(is)
        if (spcore(ist,is)) then
          evalcr(ist,jas)=evalcr(ist,ias)
          rwfcr(1:nrs,1:2,ist,jas)=rwfcr(1:nrs,1:2,ist,ias)
        end if
      end do
      rhocr(1:nr,jas,1:nspncr)=rhocr(1:nr,ias,1:nspncr)
    end if
  end do
! end loop over species and atoms
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine
!EOC

