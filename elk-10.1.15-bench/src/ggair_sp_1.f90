
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_sp_1
! !INTERFACE:
subroutine ggair_sp_1(rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
! !INPUT/OUTPUT PARAMETERS:
!   rhoup : spin-up density (in,real(ngtot))
!   rhodn : spin-down density (in,real(ngtot))
!   grho  : |grad rho| (out,real(ngtot))
!   gup   : |grad rhoup| (out,real(ngtot))
!   gdn   : |grad rhodn| (out,real(ngtot))
!   g2up  : grad^2 rhoup (out,real(ngtot))
!   g2dn  : grad^2 rhodn (out,real(ngtot))
!   g3rho : (grad rho).(grad |grad rho|) (out,real(ngtot))
!   g3up  : (grad rhoup).(grad |grad rhoup|) (out,real(ngtot))
!   g3dn  : (grad rhodn).(grad |grad rhodn|) (out,real(ngtot))
! !DESCRIPTION:
!   Computes $|\nabla\rho|$, $|\nabla\rho^{\uparrow}|$,
!   $|\nabla\rho^{\downarrow}|$, $\nabla^2\rho^{\uparrow}$,
!   $\nabla^2\rho^{\downarrow}$, $\nabla\rho\cdot(\nabla|\nabla\rho|)$,
!   $\nabla\rho^{\uparrow}\cdot(\nabla|\nabla\rho^{\uparrow}|)$ and
!   $\nabla\rho^{\downarrow}\cdot(\nabla|\nabla\rho^{\downarrow}|)$ for the
!   interstitial charge density, as required by the generalised gradient
!   approximation functionals of type 1 for spin-polarised densities. See
!   routines {\tt potxc} and {\tt modxcifc}.
!
! !REVISION HISTORY:
!   Created October 2004 (JKD)
!   Simplified and improved, October 2009 (JKD)
!EOP
!BOC
use modmain
implicit none
! arguments
real(8), intent(in) :: rhoup(ngtot),rhodn(ngtot)
real(8), intent(out) :: grho(ngtot),gup(ngtot),gdn(ngtot)
real(8), intent(out) :: g2up(ngtot),g2dn(ngtot)
real(8), intent(out) :: g3rho(ngtot),g3up(ngtot),g3dn(ngtot)
! local variables
integer ig,ifg,i
! allocatable arrays
real(8), allocatable :: gvup(:,:),gvdn(:,:),rfir(:)
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(gvup(ngtot,3),gvdn(ngtot,3),rfir(ngtot))
allocate(zfft1(nfgrz),zfft2(nfgrz))
!----------------!
!     rho up     !
!----------------!
call rzfftifc(3,ngridg,-1,rhoup,zfft1)
! |grad rhoup|
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    if (ig <= ngvc) then
      zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
    else
      zfft2(ifg)=0.d0
    end if
  end do
  call rzfftifc(3,ngridg,1,gvup(:,i),zfft2)
end do
gup(:)=sqrt(gvup(:,1)**2+gvup(:,2)**2+gvup(:,3)**2)
! grad^2 rhoup
do ifg=1,nfgrz
  ig=igrzf(ifg)
  if (ig <= ngvc) then
    zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
  else
    zfft2(ifg)=0.d0
  end if
end do
call rzfftifc(3,ngridg,1,g2up,zfft2)
! (grad rhoup).(grad |grad rhoup|)
call rzfftifc(3,ngridg,-1,gup,zfft1)
g3up(:)=0.d0
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    if (ig <= ngvec) then
      zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
    else
      zfft2(ifg)=0.d0
    end if
  end do
  call rzfftifc(3,ngridg,1,rfir,zfft2)
  g3up(:)=g3up(:)+gvup(:,i)*rfir(:)
end do
!------------------!
!     rho down     !
!------------------!
call rzfftifc(3,ngridg,-1,rhodn,zfft1)
! |grad rhodn|
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    if (ig <= ngvc) then
      zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
    else
      zfft2(ifg)=0.d0
    end if
  end do
  call rzfftifc(3,ngridg,1,gvdn(:,i),zfft2)
end do
gdn(:)=sqrt(gvdn(:,1)**2+gvdn(:,2)**2+gvdn(:,3)**2)
! grad^2 rhodn
do ifg=1,nfgrz
  ig=igrzf(ifg)
  if (ig <= ngvc) then
    zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
  else
    zfft2(ifg)=0.d0
  end if
end do
call rzfftifc(3,ngridg,1,g2dn,zfft2)
! (grad rhodn).(grad |grad rhodn|)
call rzfftifc(3,ngridg,-1,gdn,zfft1)
g3dn(:)=0.d0
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    if (ig <= ngvec) then
      zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
    else
      zfft2(ifg)=0.d0
    end if
  end do
  call rzfftifc(3,ngridg,1,rfir,zfft2)
  g3dn(:)=g3dn(:)+gvdn(:,i)*rfir(:)
end do
!-------------!
!     rho     !
!-------------!
! |grad rho|
grho(:)=sqrt((gvup(:,1)+gvdn(:,1))**2 &
            +(gvup(:,2)+gvdn(:,2))**2 &
            +(gvup(:,3)+gvdn(:,3))**2)
! (grad rho).(grad |grad rho|)
call rzfftifc(3,ngridg,-1,grho,zfft1)
g3rho(:)=0.d0
do i=1,3
  do ifg=1,nfgrz
    ig=igrzf(ifg)
    if (ig <= ngvec) then
      zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
    else
      zfft2(ifg)=0.d0
    end if
  end do
  call rzfftifc(3,ngridg,1,rfir,zfft2)
  g3rho(:)=g3rho(:)+(gvup(:,i)+gvdn(:,i))*rfir(:)
end do
deallocate(gvup,gvdn,rfir,zfft1,zfft2)
end subroutine
!EOC

