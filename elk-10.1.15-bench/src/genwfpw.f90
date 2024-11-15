
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genwfpw(vpl,ngp,igpig,vgpl,vgpc,gpc,sfacgp,nhp,vhpc,hpc,sfachp,wfpw)
use modmain
use modpw
implicit none
! arguments
real(8), intent(in) :: vpl(3)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
real(8), intent(in) :: vgpl(3,ngkmax,nspnfv),vgpc(3,ngkmax,nspnfv)
real(8), intent(in) :: gpc(ngkmax,nspnfv)
complex(8), intent(in) :: sfacgp(ngkmax,natmtot,nspnfv)
integer, intent(in) :: nhp(nspnfv)
real(8), intent(in) :: vhpc(3,nhkmax,nspnfv),hpc(nhkmax,nspnfv)
complex(8), intent(in) :: sfachp(nhkmax,natmtot,nspnfv)
complex(8), intent(out) :: wfpw(nhkmax,nspinor,nstsv)
! local variables
integer ispn0,ispn1,ispn,jspn
integer ist,is,ia,ias,igp,ihp
integer nrc,nrci,irco,irc
integer l,lm,n,i
real(8) v1,v2,v3,t0,t1
complex(8) zsm,z1,z2,z3
! automatic arrays
real(8) jl(0:lmaxo,nrcmtmax)
complex(8) ylm(lmmaxo)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgp(:,:,:)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv))
allocate(wfgp(ngkmax,nspinor,nstsv))
! get the eigenvectors from file
call getevecfv(filext,0,vpl,vgpl,evecfv)
call getevecsv(filext,0,vpl,evecsv)
! find the matching coefficients
do ispn=1,nspnfv
  call match(ngp(ispn),vgpc(:,:,ispn),gpc(:,ispn),sfacgp(:,:,ispn), &
   apwalm(:,:,:,:,ispn))
end do
! calculate the second-variational wavefunctions for all states
call genwfsv(.true.,.true.,nstsv,[0],ngridg,igfft,ngp,igpig,apwalm,evecfv, &
 evecsv,wfmt,ngkmax,wfgp)
deallocate(apwalm,evecfv,evecsv)
! zero the plane wave coefficients
wfpw(:,:,:)=0.d0
!---------------------------!
!     interstitial part     !
!---------------------------!
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
  i=1
  do ihp=1,nhp(jspn)
    v1=vhpc(1,ihp,jspn); v2=vhpc(2,ihp,jspn); v3=vhpc(3,ihp,jspn)
    do igp=i,ngp(jspn)
      t1=abs(v1-vgpc(1,igp,jspn)) &
        +abs(v2-vgpc(2,igp,jspn)) &
        +abs(v3-vgpc(3,igp,jspn))
      if (t1 < epslat) then
        wfpw(ihp,ispn0:ispn1,:)=wfgp(igp,ispn0:ispn1,:)
        if (igp == i) i=i+1
        exit
      end if
    end do
  end do
end do
!-------------------------!
!     muffin-tin part     !
!-------------------------!
t0=1.d0/sqrt(omega)
! remove continuation of interstitial function into muffin-tin
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
! loop over G+p-vectors
  do igp=1,ngp(jspn)
! generate the conjugate spherical harmonics 4π iˡ Yₗₘ(G+p)*
    call genylmv(.true.,lmaxo,vgpc(:,igp,jspn),ylm)
    ylm(:)=conjg(ylm(:))
    do is=1,nspecies
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      irco=nrci+1
! generate spherical Bessel functions
      do irc=1,nrci
        t1=gpc(igp,jspn)*rcmt(irc,is)
        call sbessel(lmaxi,t1,jl(:,irc))
      end do
      do irc=irco,nrc
        t1=gpc(igp,jspn)*rcmt(irc,is)
        call sbessel(lmaxo,t1,jl(:,irc))
      end do
! loop over atoms
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        z1=t0*sfacgp(igp,ias,jspn)
        do ist=1,nstsv
          do ispn=ispn0,ispn1
            z2=z1*wfgp(igp,ispn,ist)
            i=1
            do irc=1,nrci
              do l=0,lmaxi
                n=2*l
                lm=l**2+1
                z3=jl(l,irc)*z2
                wfmt(i:i+n,ias,ispn,ist)=wfmt(i:i+n,ias,ispn,ist) &
                 -z3*ylm(lm:lm+n)
                i=i+n+1
              end do
            end do
            do irc=irco,nrc
              do l=0,lmaxo
                n=2*l
                lm=l**2+1
                z3=jl(l,irc)*z2
                wfmt(i:i+n,ias,ispn,ist)=wfmt(i:i+n,ias,ispn,ist) &
                 -z3*ylm(lm:lm+n)
                i=i+n+1
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do
! Fourier transform the muffin-tin wavefunctions
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
! loop over H+p-vectors
  do ihp=1,nhp(jspn)
! generate the spherical harmonics 4π (-i)ˡ Yₗₘ(H+p)
    call genylmv(.true.,lmaxo,vhpc(:,ihp,jspn),ylm)
    do is=1,nspecies
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      irco=nrci+1
! generate spherical Bessel functions
      do irc=1,nrci
        t1=hpc(ihp,jspn)*rcmt(irc,is)
        call sbessel(lmaxi,t1,jl(:,irc))
      end do
      do irc=irco,nrc
        t1=hpc(ihp,jspn)*rcmt(irc,is)
        call sbessel(lmaxo,t1,jl(:,irc))
      end do
      do ia=1,natoms(is)
        ias=idxas(ia,is)
! conjugate structure factor
        z1=t0*conjg(sfachp(ihp,ias,jspn))
! loop over states
        do ist=1,nstsv
          do ispn=ispn0,ispn1
            zsm=0.d0
            i=1
            do irc=1,nrci
              z2=jl(0,irc)*wfmt(i,ias,ispn,ist)*ylm(1)
              i=i+1
              do l=1,lmaxi
                n=2*l
                lm=l**2+1
                z2=z2+jl(l,irc)*sum(wfmt(i:i+n,ias,ispn,ist)*ylm(lm:lm+n))
                i=i+n+1
              end do
              zsm=zsm+wrcmt(irc,is)*z2
            end do
            do irc=irco,nrc
              z2=jl(0,irc)*wfmt(i,ias,ispn,ist)*ylm(1)
              i=i+1
              do l=1,lmaxo
                n=2*l
                lm=l**2+1
                z2=z2+jl(l,irc)*sum(wfmt(i:i+n,ias,ispn,ist)*ylm(lm:lm+n))
                i=i+n+1
              end do
              zsm=zsm+wrcmt(irc,is)*z2
            end do
! add to the H+p wavefunction
            wfpw(ihp,ispn,ist)=wfpw(ihp,ispn,ist)+z1*zsm
          end do
        end do
      end do
    end do
  end do
end do
deallocate(wfmt,wfgp)
end subroutine

