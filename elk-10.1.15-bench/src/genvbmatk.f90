
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvbmatk(vmt,vir,bmt,bir,ngp,igpig,wfmt,ld,wfgp,vbmat)
use modmain
use moddftu
use modomp
implicit none
! arguments
! the potential and field are multiplied by the radial integration weights in
! the muffin-tin and by the characteristic function in the interstitial region
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtc)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtc,ndmag)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(4), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
integer, intent(in) :: ld
complex(4), intent(in) :: wfgp(ld,nspinor,nstsv)
complex(8), intent(out) :: vbmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn
integer is,ias,nrc,nrci,nrco
integer npc,npc2,ipco
integer n,igp,nthd
! automatic arrays
complex(4) wfmt1(npcmtmax,2),wfmt2(npcmtmax,2)
complex(4) wfir1(ngtc,nspinor),wfir2(ngtc,nspinor),c(ngkmax)
! external functions
real(4), external :: sdot
complex(4), external :: cdotc
! zero the upper triangular matrix elements
do jst=1,nstsv
  vbmat(1:jst,jst)=0.d0
end do
call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,wfmt2,wfir1,wfir2,c) &
!$OMP PRIVATE(ias,is,nrc,nrci,nrco) &
!$OMP PRIVATE(npc,npc2,ipco,ist,jst) &
!$OMP PRIVATE(ispn,jspn,igp,n) &
!$OMP NUM_THREADS(nthd)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  nrco=nrc-nrci
  npc=npcmt(is)
  npc2=npc*2
  ipco=npcmti(is)+1
!$OMP DO
  do jst=1,nstsv
! apply local potential and magnetic field to spinor wavefunction
    if (ncmag) then
! non-collinear case
      call vbmk1(npc,vmt(:,ias),bmt(:,ias,1),bmt(:,ias,2),bmt(:,ias,3), &
       wfmt(:,ias,1,jst),wfmt(:,ias,2,jst),wfmt1,wfmt1(:,2))
    else
! collinear case
      call vbmk2(npc,vmt(:,ias),bmt(:,ias,1),wfmt(:,ias,1,jst), &
       wfmt(:,ias,2,jst),wfmt1,wfmt1(:,2))
    end if
! apply muffin-tin DFT+U potential matrix if required
    if (dftu /= 0) then
      if (any(tvmmt(0:lmaxdm,ias))) then
! multiply wavefunction by radial integration weights
        wfmt2(1:npc,1)=wfmt(1:npc,ias,1,jst)
        wfmt2(1:npc,2)=wfmt(1:npc,ias,2,jst)
        call cfcmtwr(nrc,nrci,wrcmt(:,is),wfmt2(1,1))
        call cfcmtwr(nrc,nrci,wrcmt(:,is),wfmt2(1,2))
        call cgemm('N','N',lmmaxi,nrci,lmmaxi,cone,vmatmti(1,1,1,1,ias), &
         lmmaxi,wfmt2(1,1),lmmaxi,cone,wfmt1(1,1),lmmaxi)
        call cgemm('N','N',lmmaxo,nrco,lmmaxo,cone,vmatmto(1,1,1,1,ias), &
         lmmaxo,wfmt2(ipco,1),lmmaxo,cone,wfmt1(ipco,1),lmmaxo)
        call cgemm('N','N',lmmaxi,nrci,lmmaxi,cone,vmatmti(1,1,2,2,ias), &
         lmmaxi,wfmt2(1,2),lmmaxi,cone,wfmt1(1,2),lmmaxi)
        call cgemm('N','N',lmmaxo,nrco,lmmaxo,cone,vmatmto(1,1,2,2,ias), &
         lmmaxo,wfmt2(ipco,2),lmmaxo,cone,wfmt1(ipco,2),lmmaxo)
        if (ncmag) then
          call cgemm('N','N',lmmaxi,nrci,lmmaxi,cone,vmatmti(1,1,1,2,ias), &
           lmmaxi,wfmt2(1,2),lmmaxi,cone,wfmt1(1,1),lmmaxi)
          call cgemm('N','N',lmmaxo,nrco,lmmaxo,cone,vmatmto(1,1,1,2,ias), &
           lmmaxo,wfmt2(ipco,2),lmmaxo,cone,wfmt1(ipco,1),lmmaxo)
          call cgemm('N','N',lmmaxi,nrci,lmmaxi,cone,vmatmti(1,1,2,1,ias), &
           lmmaxi,wfmt2(1,1),lmmaxi,cone,wfmt1(1,2),lmmaxi)
          call cgemm('N','N',lmmaxo,nrco,lmmaxo,cone,vmatmto(1,1,2,1,ias), &
           lmmaxo,wfmt2(ipco,1),lmmaxo,cone,wfmt1(ipco,2),lmmaxo)
        end if
      end if
    end if
! compute the inner products
    do ist=1,jst-1
      vbmat(ist,jst)=vbmat(ist,jst) &
       +cdotc(npc,wfmt(1,ias,1,ist),1,wfmt1(1,1),1) &
       +cdotc(npc,wfmt(1,ias,2,ist),1,wfmt1(1,2),1)
    end do
    vbmat(jst,jst)=vbmat(jst,jst) &
     +sdot(npc2,wfmt(1,ias,1,jst),1,wfmt1(1,1),1) &
     +sdot(npc2,wfmt(1,ias,2,jst),1,wfmt1(1,2),1)
  end do
!$OMP END DO
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
!$OMP DO
do jst=1,nstsv
! Fourier transform wavefunction to real-space
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    wfir1(1:ngtc,ispn)=0.e0
    do igp=1,ngp(jspn)
      wfir1(igfc(igpig(igp,jspn)),ispn)=wfgp(igp,ispn,jst)
    end do
    call cfftifc(3,ngdgc,1,wfir1(:,ispn))
  end do
! apply local potential and magnetic field to spinor wavefunction
  if (ncmag) then
! non-collinear case
    call vbmk1(ngtc,vir,bir,bir(:,2),bir(:,3),wfir1,wfir1(:,2),wfir2,wfir2(:,2))
  else
! collinear case
    call vbmk2(ngtc,vir,bir,wfir1,wfir1(:,2),wfir2,wfir2(:,2))
  end if
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    n=ngp(jspn)
! Fourier transform to G+p-space
    call cfftifc(3,ngdgc,-1,wfir2(:,ispn))
    do igp=1,n
      c(igp)=wfir2(igfc(igpig(igp,jspn)),ispn)
    end do
    do ist=1,jst-1
      vbmat(ist,jst)=vbmat(ist,jst)+cdotc(n,wfgp(1,ispn,ist),1,c,1)
    end do
    vbmat(jst,jst)=vbmat(jst,jst)+sdot(2*n,wfgp(1,ispn,jst),1,c,1)
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vbmat(ist,jst)=conjg(vbmat(jst,ist))
  end do
end do
return

contains

pure subroutine vbmk1(n,v,b1,b2,b3,wf11,wf12,wf21,wf22)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: v(n),b1(n),b2(n),b3(n)
complex(4), intent(in) :: wf11(n),wf12(n)
complex(4), intent(out) :: wf21(n),wf22(n)
! local variables
integer i
!$OMP SIMD SIMDLEN(8)
do i=1,n
  wf21(i)=(v(i)+b3(i))*wf11(i)+cmplx(b1(i),-b2(i),8)*wf12(i)
  wf22(i)=(v(i)-b3(i))*wf12(i)+cmplx(b1(i),b2(i),8)*wf11(i)
end do
end subroutine

pure subroutine vbmk2(n,v,b,wf11,wf12,wf21,wf22)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: v(n),b(n)
complex(4), intent(in) :: wf11(n),wf12(n)
complex(4), intent(out) :: wf21(n),wf22(n)
! local variables
integer i
!$OMP SIMD SIMDLEN(8)
do i=1,n
  wf21(i)=(v(i)+b(i))*wf11(i)
  wf22(i)=(v(i)-b(i))*wf12(i)
end do
end subroutine

end subroutine

