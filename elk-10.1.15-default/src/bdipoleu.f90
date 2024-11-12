
! Copyright (C) 2024 Wenhan Chen, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bdipoleu
use modmain
use modulr
implicit none
! local variables
integer idm,iq,ifq
real(8) cb,t1
complex(8) zv(3),z1
! automatic arrays
real(8) rfft(nqpt)
complex(8) zvq(nfqrz,3)
if (.not.ncmag) then
  write(*,*)
  write(*,'("Error(bdipoleu): non-collinear magnetism required for inclusion &
   &of the dipole field")')
  write(*,*)
  stop
end if
! prefactor for the electron spin dipole magnetic field
cb=gfacte/(4.d0*solsc)
! Fourier transform the R-dependent magnetisation to Q-space
do idm=1,3
  rfft(1:nqpt)=momtotru(idm,1:nqpt)/omega
  call rzfftifc(3,ngridq,-1,rfft,zvq(:,idm))
end do
do ifq=1,nfqrz
  iq=iqrzf(ifq)
! compute B(Q) = 4π (gₛ/4c)² [ m(Q) - Q (Q⋅m(Q))/Q² ]
  zv(1:3)=zvq(ifq,1:3)
  if (iq > 1) then
    z1=vqc(1,iq)*zv(1)+vqc(2,iq)*zv(2)+vqc(3,iq)*zv(3)
    t1=vqc(1,iq)**2+vqc(2,iq)**2+vqc(3,iq)**2
    z1=z1/t1
    zv(1:3)=zv(1:3)-z1*vqc(:,iq)
  end if
  t1=bdipscf*fourpi*cb**2
  bdipq(1:3,ifq)=t1*zv(1:3)
end do
end subroutine

