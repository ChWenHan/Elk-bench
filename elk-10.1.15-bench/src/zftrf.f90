
! Copyright (C) 2010 Alexey I. Baranov.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zftrf
! !INTERFACE:
subroutine zftrf(npv,ivp,vpc,rfmt,rfir,zfp)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   npv  : number of P-vectors (in,integer)
!   ivp  : integer coordinates of the P-vectors (in,integer(3,npv))
!   vpc  : P-vectors in Cartesian coordinates (in,real(3,npv))
!   rfmt : real muffin-tin function (in,real(npmtmax,natmtot))
!   rfir : real interstitial function (in,real(ngtot))
!   zfp  : Fourier expansion coefficients of the real-space function
!          (out,complex(npv))
! !DESCRIPTION:
!   Given a real function periodic in the unit cell, $f({\bf r})$, this routine
!   calculates its complex Fourier expansion coefficients:
!   $$ f({\bf P})=\frac{1}{\Omega}\int d^3r\,f({\bf r})\tilde{\Theta}({\bf r})
!    e^{-i{\bf P}\cdot{\bf r}}
!    +\frac{4\pi}{\Omega}\sum_{\alpha}e^{-i{\bf P}\cdot{\bf R}_{\alpha}}
!    \sum_{lm}(-i)^l Y_{lm}(\hat{\bf P})
!    \int_{0}^{R_{\alpha}}dr\,r^2 j_{l}(|{\bf P}|r)f_{lm}^{\alpha}(r), $$
!   where $\tilde{\Theta}$ is the smooth characteristic function of the
!   interstitial region, $\Omega$ is the unit cell volume and $R_{\alpha}$ is
!   the muffin-tin radius of atom $\alpha$.
!
! !REVISION HISTORY:
!   Created July 2010 (Alexey I. Baranov)
!   Modified, November 2010 (JKD)
!   Optimised, May 2024 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: npv,ivp(3,npv)
real(8), intent(in) :: vpc(3,npv)
real(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
complex(8), intent(out) :: zfp(npv)
! local variables
integer is,ia,ias,ip,ig,ifg
integer nrc,nrci,irco,irc
integer l,lm,n,i
real(8) p,t0,t1
complex(8) zsm,z1,z2
! automatic arrays
real(8) jl(0:lmaxo,nrcmtmax),rfmt1(npcmtmax)
complex(8) ylm(lmmaxo)
! allocatable arrays
complex(8), allocatable :: zfft(:),zfmt(:,:)
allocate(zfft(ngtot),zfmt(npcmtmax,natmtot))
! zero the coefficients
zfp(1:npv)=0.d0
!---------------------------!
!     interstitial part     !
!---------------------------!
! Fourier transform to G-space
zfft(1:ngtot)=rfir(1:ngtot)
call zfftifc(3,ngridg,-1,zfft)
! find coefficients for all required input vectors
do ip=1,npv
  if ((ivp(1,ip) >= intgv(1,1)).and.(ivp(1,ip) <= intgv(2,1)).and. &
      (ivp(2,ip) >= intgv(1,2)).and.(ivp(2,ip) <= intgv(2,2)).and. &
      (ivp(3,ip) >= intgv(1,3)).and.(ivp(3,ip) <= intgv(2,3))) then
    ig=ivgig(ivp(1,ip),ivp(2,ip),ivp(3,ip))
    zfp(ip)=zfft(igfft(ig))
  end if
end do
!-------------------------!
!     muffin-tin part     !
!-------------------------!
! convert function from real to complex spherical harmonic expansion on coarse
! radial mesh
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  call rfmtftoc(nrc,nrci,rfmt(:,ias),rfmt1)
  call rtozfmt(nrc,nrci,rfmt1,zfmt(:,ias))
end do
! remove continuation of interstitial function into muffin-tin
do ig=1,ngvec
  ifg=igfft(ig)
! conjugate spherical harmonics 4π iˡ Yₗₘ(G)*
  call genylmv(.true.,lmaxo,vgc(:,ig),ylm)
  ylm(1:lmmaxo)=conjg(ylm(1:lmmaxo))
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    irco=nrci+1
! generate spherical Bessel functions
    do irc=1,nrci
      t1=gc(ig)*rcmt(irc,is)
      call sbessel(lmaxi,t1,jl(:,irc))
    end do
    do irc=irco,nrc
      t1=gc(ig)*rcmt(irc,is)
      call sbessel(lmaxo,t1,jl(:,irc))
    end do
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! structure factor
      t1=vgc(1,ig)*atposc(1,ia,is) &
        +vgc(2,ig)*atposc(2,ia,is) &
        +vgc(3,ig)*atposc(3,ia,is)
      z1=zfft(ifg)*cmplx(cos(t1),sin(t1),8)
      i=1
      do irc=1,nrci
        do l=0,lmaxi
          n=2*l
          lm=l**2+1
          z2=z1*jl(l,irc)
          zfmt(i:i+n,ias)=zfmt(i:i+n,ias)-z2*ylm(lm:lm+n)
          i=i+n+1
        end do
      end do
      do irc=irco,nrc
        do l=0,lmaxo
          n=2*l
          lm=l**2+1
          z2=z1*jl(l,irc)
          zfmt(i:i+n,ias)=zfmt(i:i+n,ias)-z2*ylm(lm:lm+n)
          i=i+n+1
        end do
      end do
    end do
  end do
end do
t0=1.d0/omega
! loop over input P-vectors
do ip=1,npv
! length of P-vector
  p=sqrt(vpc(1,ip)**2+vpc(2,ip)**2+vpc(3,ip)**2)
! spherical harmonics 4π (-i)ˡ Yₗₘ(P)
  call genylmv(.true.,lmaxo,vpc(:,ip),ylm)
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    irco=nrci+1
! generate spherical Bessel functions
    do irc=1,nrci
      t1=p*rcmt(irc,is)
      call sbessel(lmaxi,t1,jl(:,irc))
    end do
    do irc=irco,nrc
      t1=p*rcmt(irc,is)
      call sbessel(lmaxo,t1,jl(:,irc))
    end do
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      zsm=0.d0
      i=1
      do irc=1,nrci
        z1=jl(0,irc)*zfmt(i,ias)*ylm(1)
        i=i+1
        do l=1,lmaxi
          n=2*l
          lm=l**2+1
          z1=z1+jl(l,irc)*sum(zfmt(i:i+n,ias)*ylm(lm:lm+n))
          i=i+n+1
        end do
        zsm=zsm+wrcmt(irc,is)*z1
      end do
      do irc=irco,nrc
        z1=jl(0,irc)*zfmt(i,ias)*ylm(1)
        i=i+1
        do l=1,lmaxo
          n=2*l
          lm=l**2+1
          z1=z1+jl(l,irc)*sum(zfmt(i:i+n,ias)*ylm(lm:lm+n))
          i=i+n+1
        end do
        zsm=zsm+wrcmt(irc,is)*z1
      end do
! conjugate structure factor
      t1=vpc(1,ip)*atposc(1,ia,is) &
        +vpc(2,ip)*atposc(2,ia,is) &
        +vpc(3,ip)*atposc(3,ia,is)
      z1=t0*cmplx(cos(t1),-sin(t1),8)
      zfp(ip)=zfp(ip)+z1*zsm
    end do
  end do
end do
deallocate(zfft,zfmt)
end subroutine
! EOC

