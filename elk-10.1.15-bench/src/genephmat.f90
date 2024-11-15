
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genephmat(iq,ik,de,a,dvmt,dvir,ephmat)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq,ik
real(8), intent(in) :: de
complex(8), intent(in) :: a(nbph,nbph)
complex(8), intent(in) :: dvmt(npcmtmax,natmtot,nbph),dvir(ngtc,nbph)
complex(8), intent(out) :: ephmat(nstsv,nstsv,nbph)
! local variables
integer jk,jkq,isym,ld
integer nst,nstq,ist,jst
integer ispn,jspn,is,ias
integer npc,nq,igp,i,j,l
real(8) vpql(3)
! automatic arrays
integer idx(nstsv),idxq(nstsv)
integer ngp(nspnfv),ngpq(nspnfv)
complex(4) cfmt1(npcmtmax),cfmt2(npcmtmax),c(ngkmax)
complex(8) x(nbph)
! allocatable arrays
integer, allocatable :: igpig(:,:),igpqig(:,:)
complex(4), allocatable :: wfmt(:,:,:,:),wfgp(:,:,:)
complex(4), allocatable :: wfmtq(:,:,:,:),wfgpq(:,:,:)
complex(4), allocatable :: wfir1(:),wfir2(:)
! equivalent reduced k-point
jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! k+q-vector in lattice coordinates
vpql(1:3)=vkl(1:3,ik)+vql(1:3,iq)
! find reduced k-point index corresponding to k+q
call findkpt(vpql,isym,jkq)
! index to states in energy window around Fermi energy
nst=0
do ist=1,nstsv
  if (abs(evalsv(ist,jk)-efermi) > de) cycle
  nst=nst+1
  idx(nst)=ist
end do
nstq=0
do ist=1,nstsv
  if (abs(evalsv(ist,jkq)-efermi) > de) cycle
  nstq=nstq+1
  idxq(nstq)=ist
end do
! generate the second-variational wavefunctions for all states at k and k+q
allocate(igpig(ngkmax,nspnfv))
allocate(wfmt(npcmtmax,natmtot,nspinor,nst),wfgp(ngkmax,nspinor,nst))
call genwfsvp_sp(.false.,.true.,nst,idx,ngridg,igfft,vkl(:,ik),ngp,igpig,wfmt, &
 ngkmax,wfgp)
allocate(igpqig(ngkmax,nspnfv))
allocate(wfmtq(npcmtmax,natmtot,nspinor,nstq),wfgpq(ngkmax,nspinor,nstq))
call genwfsvp_sp(.false.,.true.,nstq,idxq,ngridg,igfft,vpql,ngpq,igpqig,wfmtq, &
 ngkmax,wfgpq)
! zero the electron-phonon coupling matrix elements
ephmat(1:nstsv,1:nstsv,1:nbph)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do j=1,nst
  jst=idx(j)
  do i=1,nstq
    ist=idxq(i)
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      if (spinpol) then
        cfmt1(1:npc)=wfmtq(1:npc,ias,1,i)*conjg(wfmt(1:npc,ias,1,j)) &
                    +wfmtq(1:npc,ias,2,i)*conjg(wfmt(1:npc,ias,2,j))
      else
        cfmt1(1:npc)=wfmtq(1:npc,ias,1,i)*conjg(wfmt(1:npc,ias,1,j))
      end if
      call cfsht(nrcmt(is),nrcmti(is),cfmt1,cfmt2)
      do l=1,nbph
        ephmat(ist,jst,l)=ephmat(ist,jst,l) &
         +dot_product(cfmt2(1:npc),dvmt(1:npc,ias,l))
      end do
    end do
  end do
end do
deallocate(wfmt,wfmtq)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(wfir1(ngtc),wfir2(ngtc))
do j=1,nst
  jst=idx(j)
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    nq=ngpq(jspn)
! Fourier transform wavefunction to real-space
    wfir1(1:ngtc)=0.e0
    do igp=1,ngp(jspn)
      wfir1(igfc(igpig(igp,jspn)))=wfgp(igp,ispn,j)
    end do
    call cfftifc(3,ngdgc,1,wfir1)
    do l=1,nbph
! apply potential derivative to wavefunction
      wfir2(1:ngtc)=dvir(1:ngtc,l)*wfir1(1:ngtc)
! Fourier transform to G+p+q-space
      call cfftifc(3,ngdgc,-1,wfir2)
      do igp=1,nq
        c(igp)=wfir2(igfc(igpqig(igp,jspn)))
      end do
      do i=1,nstq
        ist=idxq(i)
! compute inner product
        ephmat(ist,jst,l)=ephmat(ist,jst,l) &
         +dot_product(wfgpq(1:nq,ispn,i),c(1:nq))
      end do
    end do
  end do
end do
deallocate(wfir1,wfir2)
! convert to phonon coordinates
ld=nstsv**2
do i=1,nstq
  ist=idxq(i)
  do j=1,nst
    jst=idx(j)
    x(1:nbph)=ephmat(ist,jst,1:nbph)
    call zgemv('T',nbph,nbph,zone,a,nbph,x,1,zzero,ephmat(ist,jst,1),ld)
  end do
end do
deallocate(igpig,igpqig,wfgp,wfgpq)
end subroutine

