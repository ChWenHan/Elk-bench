
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdisp
use modmain
use modphonon
implicit none
! local variables
integer iq,i,iv
real(8) wmin,wmax,t1
! allocatable arrays
real(8), allocatable :: wq(:,:)
complex(8), allocatable :: dq(:,:),ev(:,:)
! initialise universal variables
call init0
call init2
call initph
! allocate local arrays
allocate(wq(nbph,npp1d),dq(nbph,nbph),ev(nbph,nbph))
! generate a set of q-point vectors along a path in the Brillouin zone
call plotpt1d(bvec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
! compute the phonon frequencies along the path
do iq=ip01d,npp1d
! compute the dynamical matrix at this particular q-point
  call dynrtoq(vplp1d(:,iq),dynr,dq)
! find the phonon frequencies and eigenvectors
  call dynev(dq,wq(:,iq),ev)
end do
! find the minimum and maximum phonon frequencies
wmin=minval(wq(1,:))
wmax=maxval(wq(nbph,:))
t1=(wmax-wmin)*0.5d0
wmin=wmin-t1
wmax=wmax+t1
! output the vertex location lines
open(50,file='PHDLINES.OUT',form='FORMATTED')
do iv=1,nvp1d
  write(50,'(2G18.10)') dvp1d(iv),wmin
  write(50,'(2G18.10)') dvp1d(iv),wmax
  write(50,*)
end do
close(50)
! output the phonon dispersion
open(50,file='PHDISP.OUT',form='FORMATTED')
do i=1,nbph
  do iq=ip01d,npp1d
    write(50,'(2G18.10)') dpp1d(iq),wq(i,iq)
  end do
  write(50,*)
end do
close(50)
write(*,*)
write(*,'("Info(phdisp):")')
write(*,'(" phonon dispersion written to PHDISP.OUT")')
write(*,'(" vertex location lines written to PHDLINES.OUT")')
deallocate(dq,ev)
end subroutine

