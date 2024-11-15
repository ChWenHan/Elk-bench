
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl,
! F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnsv(ngp,igpig,vgpc,apwalm,evalfv,evecfv,evalsvp,evecsv)
use modmain
use moddftu
use modomp
implicit none
! arguments
integer, intent(in) :: ngp,igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(in) :: evalfv(nstfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
real(8), intent(out) :: evalsvp(nstsv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
logical todsb,socz
integer ld,ist,jst,ispn,is,ias
integer nrc,nrci,nrco,irco,irc
integer l,lm,nm,npc,npc2,npci,ipco
integer ngp2,igp,i0,i,j,nthd
real(8) ca,a(3),asp(3,3),t1
real(8) ts0,ts1
complex(8) z1,z2,z3
complex(4) c1
! automatic arrays
complex(8) wfmt2(npcmtmax),wfmt4(npcmtmax,3)
complex(8) wfmt31(npcmtmax),wfmt32(npcmtmax),wfmt33(npcmtmax)
complex(4) wfmt5(npcmtmax),wfgp1(ngkmax),wfgp2(ngkmax),wfgp3(ngkmax)
complex(4) wfir1(ngtc),wfir2(ngtc)
! allocatable arrays
complex(4), allocatable :: wfmt0(:,:),wfgp0(:,:)
complex(8), allocatable :: wfmt1(:,:)
! external functions
real(4), external :: sdot
complex(4), external :: cdotc
! no calculation of second-variational eigenvectors
if (.not.tevecsv) then
  evalsvp(1:nstsv)=evalfv(1:nstsv)
  evecsv(1:nstsv,1:nstsv)=0.d0
  do i=1,nstsv
    evecsv(i,i)=1.d0
  end do
  return
end if
call timesec(ts0)
! coupling constant of the external A-field (-1/c)
ca=-1.d0/solsc
if (tafield) a(1:3)=ca*afieldc(1:3)
if (tafsp) asp(1:3,1:3)=ca*afspc(1:3,1:3)
! check if the off-diagonal spin block of the Hamiltonian is required
if (spinpol.and.(ncmag.or.spinorb)) then
  todsb=.true.
else
  todsb=.false.
end if
! special case of spin-orbit coupling and collinear magnetism
if (spinorb.and.cmagz) then
  socz=.true.
else
  socz=.false.
end if
ld=lmmaxdm*nspinor
call holdthd(nstfv,nthd)
! zero the second-variational Hamiltonian (stored in the eigenvector array)
evecsv(1:nstsv,1:nstsv)=0.d0
! set the diagonal elements equal to the first-variational eigenvalues
do ispn=1,nspinor
  do ist=1,nstfv
    i=nstfv*(ispn-1)+ist
    evecsv(i,i)=evalfv(ist)
  end do
end do
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(wfmt0(npcmtmax,nstfv),wfmt1(npcmtmax,nstfv))
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt2,wfmt31,wfmt32,wfmt33,wfmt4,wfmt5) &
!$OMP PRIVATE(ias,is,nrc,nrci,nrco,irco) &
!$OMP PRIVATE(npc,npc2,npci,ipco,ist,jst,irc) &
!$OMP PRIVATE(t1,i0,i,j,l,nm,lm,z1,z2,z3) &
!$OMP NUM_THREADS(nthd)
! begin loop over atoms
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  nrco=nrc-nrci
  irco=nrci+1
  npc=npcmt(is)
  npc2=npc*2
  npci=npcmti(is)
  ipco=npci+1
! compute the first-variational wavefunctions
!$OMP DO
  do ist=1,nstfv
    call wfmtfv(ias,ngp,apwalm(:,:,:,ias),evecfv(:,ist),wfmt1(:,ist))
! make a single-precision copy of the wavefunction
    wfmt0(1:npc,ist)=wfmt1(1:npc,ist)
  end do
!$OMP END DO
! begin loop over states
!$OMP DO
  do jst=1,nstfv
    if (spinpol) then
! convert wavefunction to spherical coordinates
      call zbsht(nrc,nrci,wfmt1(:,jst),wfmt2)
! apply Kohn-Sham effective magnetic field
      wfmt32(1:npc)=bsmt(1:npc,ias,ndmag)*wfmt2(1:npc)
! convert to spherical harmonics
      call zfsht(nrc,nrci,wfmt32,wfmt31)
! non-collinear magnetic field
      if (socz) then
        wfmt33(1:npc)=0.d0
      else if (ncmag) then
        wfmt32(1:npc)=cmplx(bsmt(1:npc,ias,1),-bsmt(1:npc,ias,2),8)*wfmt2(1:npc)
        call zfsht(nrc,nrci,wfmt32,wfmt33)
      end if
      wfmt32(1:npc)=-wfmt31(1:npc)
! apply spin-orbit coupling if required
      if (spinorb) then
        call lopzflmn(lmaxi,nrci,lmmaxi,wfmt1(1,jst),wfmt4,wfmt4(1,2), &
         wfmt4(1,3))
        call lopzflmn(lmaxo,nrco,lmmaxo,wfmt1(ipco,jst),wfmt4(ipco,1), &
         wfmt4(ipco,2),wfmt4(ipco,3))
! inner part of muffin-tin
        do irc=1,nrci
          t1=socfr(irc,ias)
          i0=lmmaxi*(irc-1)
          do i=i0+1,i0+lmmaxi
            wfmt31(i)=wfmt31(i)+t1*wfmt4(i,3)
            wfmt32(i)=wfmt32(i)-t1*wfmt4(i,3)
            wfmt33(i)=wfmt33(i)+t1*(wfmt4(i,1) &
             +cmplx(aimag(wfmt4(i,2)),-dble(wfmt4(i,2)),8))
          end do
        end do
! outer part of muffin-tin
        do irc=irco,nrc
          t1=socfr(irc,ias)
          i0=npci+lmmaxo*(irc-irco)
          do i=i0+1,i0+lmmaxo
            wfmt31(i)=wfmt31(i)+t1*wfmt4(i,3)
            wfmt32(i)=wfmt32(i)-t1*wfmt4(i,3)
            wfmt33(i)=wfmt33(i)+t1*(wfmt4(i,1) &
             +cmplx(aimag(wfmt4(i,2)),-dble(wfmt4(i,2)),8))
          end do
        end do
      end if
    else
      wfmt31(1:npc)=0.d0
    end if
! apply muffin-tin DFT+U potential matrix if required
    if (tvmatmt) then
      do l=0,lmaxdm
        if (tvmmt(l,ias)) then
          nm=2*l+1
          lm=l**2+1
          i=npci+lm
          if (l <= lmaxi) then
            call zgemm('N','N',nm,nrci,nm,zone,vmatmt(lm,1,lm,1,ias),ld, &
             wfmt1(lm,jst),lmmaxi,zone,wfmt31(lm),lmmaxi)
          end if
          call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,1,lm,1,ias),ld, &
           wfmt1(i,jst),lmmaxo,zone,wfmt31(i),lmmaxo)
          if (spinpol) then
            if (l <= lmaxi) then
              call zgemm('N','N',nm,nrci,nm,zone,vmatmt(lm,2,lm,2,ias),ld, &
               wfmt1(lm,jst),lmmaxi,zone,wfmt32(lm),lmmaxi)
            end if
            call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,2,lm,2,ias),ld, &
             wfmt1(i,jst),lmmaxo,zone,wfmt32(i),lmmaxo)
            if (todsb) then
              if (l <= lmaxi) then
                call zgemm('N','N',nm,nrci,nm,zone,vmatmt(lm,1,lm,2,ias),ld, &
                 wfmt1(lm,jst),lmmaxi,zone,wfmt33(lm),lmmaxi)
              end if
              call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,1,lm,2,ias),ld, &
               wfmt1(i,jst),lmmaxo,zone,wfmt33(i),lmmaxo)
            end if
          end if
        end if
      end do
    end if
! apply vector potential if required
    if (tafield.or.tafsp) then
      call gradzfmt(nrc,nrci,rlcmt(:,-1,is),wcrcmt(:,:,is),wfmt1(:,jst), &
       npcmtmax,wfmt4)
      if (tafield) then
        do i=1,npc
          z1=a(1)*wfmt4(i,1)+a(2)*wfmt4(i,2)+a(3)*wfmt4(i,3)
          z1=cmplx(aimag(z1),-dble(z1),8)
          wfmt31(i)=wfmt31(i)+z1
          if (spinpol) wfmt32(i)=wfmt32(i)+z1
        end do
      end if
! apply spin-dependent vector potential if required
      if (tafsp) then
        do i=1,npc
          z3=asp(1,3)*wfmt4(i,1)+asp(2,3)*wfmt4(i,2)+asp(3,3)*wfmt4(i,3)
          z3=cmplx(aimag(z3),-dble(z3),8)
          wfmt31(i)=wfmt31(i)+z3
          wfmt32(i)=wfmt32(i)-z3
          if (ncmag) then
            z1=asp(1,1)*wfmt4(i,1)+asp(2,1)*wfmt4(i,2)+asp(3,1)*wfmt4(i,3)
            z2=asp(1,2)*wfmt4(i,1)+asp(2,2)*wfmt4(i,2)+asp(3,2)*wfmt4(i,3)
            wfmt33(i)=wfmt33(i)+cmplx(aimag(z1),-dble(z1),8)-z2
          end if
        end do
      end if
    end if
! add to second-variational Hamiltonian matrix
! upper diagonal block
    call zcfmtwr(nrc,nrci,wrcmt(:,is),wfmt31,wfmt5)
    do ist=1,jst-1
      evecsv(ist,jst)=evecsv(ist,jst)+cdotc(npc,wfmt0(:,ist),1,wfmt5,1)
    end do
    evecsv(jst,jst)=evecsv(jst,jst)+sdot(npc2,wfmt0(:,jst),1,wfmt5,1)
    if (spinpol) then
      j=jst+nstfv
! lower diagonal block
      call zcfmtwr(nrc,nrci,wrcmt(:,is),wfmt32,wfmt5)
      do ist=1,jst-1
        i=ist+nstfv
        evecsv(i,j)=evecsv(i,j)+cdotc(npc,wfmt0(:,ist),1,wfmt5,1)
      end do
      evecsv(j,j)=evecsv(j,j)+sdot(npc2,wfmt0(:,jst),1,wfmt5,1)
! off-diagonal block
      if (todsb) then
        call zcfmtwr(nrc,nrci,wrcmt(:,is),wfmt33,wfmt5)
        do ist=1,nstfv
          evecsv(ist,j)=evecsv(ist,j)+cdotc(npc,wfmt0(:,ist),1,wfmt5,1)
        end do
      end if
    end if
! end loop over states
  end do
!$OMP END DO
! end loop over atoms
end do
!$OMP END PARALLEL
deallocate(wfmt0,wfmt1)
!---------------------------!
!     interstitial part     !
!---------------------------!
if (spinpol.or.tafield) then
  if (socz) todsb=.false.
  ngp2=ngp*2
  allocate(wfgp0(ngp,nstfv))
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,wfir2,wfgp1,wfgp2,wfgp3) &
!$OMP PRIVATE(ist,jst,igp,t1,c1,i,j) &
!$OMP NUM_THREADS(nthd)
! make single-precision copy of wavefunction
!$OMP DO
  do ist=1,nstfv
    wfgp0(1:ngp,ist)=evecfv(1:ngp,ist)
  end do
!$OMP END DO
! begin loop over states
!$OMP DO
  do jst=1,nstfv
    wfir1(1:ngtc)=0.e0
    do igp=1,ngp
      wfir1(igfc(igpig(igp)))=wfgp0(igp,jst)
    end do
! Fourier transform wavefunction to real-space
    call cfftifc(3,ngdgc,1,wfir1)
! multiply with magnetic field and transform to G-space
    if (spinpol) then
      wfir2(1:ngtc)=bsirc(1:ngtc,ndmag)*wfir1(1:ngtc)
      call cfftifc(3,ngdgc,-1,wfir2)
      do igp=1,ngp
        wfgp1(igp)=wfir2(igfc(igpig(igp)))
      end do
      wfgp2(1:ngp)=-wfgp1(1:ngp)
      if (ncmag) then
        wfir2(1:ngtc)=cmplx(bsirc(1:ngtc,1),-bsirc(1:ngtc,2),8)*wfir1(1:ngtc)
        call cfftifc(3,ngdgc,-1,wfir2)
        do igp=1,ngp
          wfgp3(igp)=wfir2(igfc(igpig(igp)))
        end do
      end if
    else
      wfgp1(1:ngp)=0.e0
    end if
! apply vector potential if required
    if (tafield) then
      wfir1(1:ngtc)=0.e0
      do igp=1,ngp
        t1=a(1)*vgpc(1,igp)+a(2)*vgpc(2,igp)+a(3)*vgpc(3,igp)
        wfir1(igfc(igpig(igp)))=t1*wfgp0(igp,jst)
      end do
      call cfftifc(3,ngdgc,1,wfir1)
      wfir1(1:ngtc)=wfir1(1:ngtc)*cfrc(1:ngtc)
      call cfftifc(3,ngdgc,-1,wfir1)
      do igp=1,ngp
        c1=wfir1(igfc(igpig(igp)))
        wfgp1(igp)=wfgp1(igp)+c1
        if (spinpol) wfgp2(igp)=wfgp2(igp)+c1
      end do
    end if
! apply spin-dependent vector potential if required
    if (tafsp) then
      do j=1,3
        if (sum(abs(asp(1:3,j))) < 1.d-8) cycle
        wfir1(1:ngtc)=0.e0
        do igp=1,ngp
          t1=asp(1,j)*vgpc(1,igp)+asp(2,j)*vgpc(2,igp)+asp(3,j)*vgpc(3,igp)
          wfir1(igfc(igpig(igp)))=t1*wfgp0(igp,jst)
        end do
        call cfftifc(3,ngdgc,1,wfir1)
        wfir1(1:ngtc)=wfir1(1:ngtc)*cfrc(1:ngtc)
        call cfftifc(3,ngdgc,-1,wfir1)
        if (j == 1) then
          do igp=1,ngp
            wfgp3(igp)=wfgp3(igp)+wfir1(igfc(igpig(igp)))
          end do
        else if (j == 2) then
          do igp=1,ngp
            c1=wfir1(igfc(igpig(igp)))
            wfgp3(igp)=wfgp3(igp)+cmplx(aimag(c1),-real(c1),4)
          end do
        else
          do igp=1,ngp
            c1=wfir1(igfc(igpig(igp)))
            wfgp1(igp)=wfgp1(igp)+c1
            wfgp2(igp)=wfgp2(igp)-c1
          end do
        end if
      end do
    end if
! add to second-variational Hamiltonian matrix
! upper diagonal block
    do ist=1,jst-1
      evecsv(ist,jst)=evecsv(ist,jst)+cdotc(ngp,wfgp0(:,ist),1,wfgp1,1)
    end do
    evecsv(jst,jst)=evecsv(jst,jst)+sdot(ngp2,wfgp0(:,jst),1,wfgp1,1)
    if (spinpol) then
      j=jst+nstfv
! lower diagonal block
      do ist=1,jst-1
        i=ist+nstfv
        evecsv(i,j)=evecsv(i,j)+cdotc(ngp,wfgp0(:,ist),1,wfgp2,1)
      end do
      evecsv(j,j)=evecsv(j,j)+sdot(ngp2,wfgp0(:,jst),1,wfgp2,1)
! off-diagonal block
      if (todsb) then
        do ist=1,nstfv
          evecsv(ist,j)=evecsv(ist,j)+cdotc(ngp,wfgp0(:,ist),1,wfgp3,1)
        end do
      end if
    end if
! end loop over states
  end do
!$OMP END DO
!$OMP END PARALLEL
  deallocate(wfgp0)
end if
call freethd(nthd)
if (spcpl.or.(.not.spinpol)) then
! spins are coupled; or spin-unpolarised: full diagonalisation
  call eveqnzh(nstsv,nstsv,evecsv,evalsvp)
else
! spins not coupled: block diagonalise H
  call eveqnzh(nstfv,nstsv,evecsv,evalsvp)
  evecsv(nstfv+1:nstsv,1:nstfv)=0.d0
  evecsv(1:nstfv,nstfv+1:nstsv)=0.d0
  i=nstfv+1
  call eveqnzh(nstfv,nstsv,evecsv(i,i),evalsvp(i))
end if
call timesec(ts1)
!$OMP ATOMIC
timesv=timesv+ts1-ts0
end subroutine

