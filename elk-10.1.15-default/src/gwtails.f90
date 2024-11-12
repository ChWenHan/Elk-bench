
! Copyright (C) 2018 A. Davydov, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gwtails
! !INTERFACE:
pure complex(8) function gwtails(ge)
! !USES:
use modmain
use modgw
! !INPUT/OUTPUT PARAMETERS:
!   ge : Green's function at the Matsubara end points (in,complex(4))
! !DESCRIPTION:
!   Sums the tails of the Green's function over the Matsubara frequencies as
!   part of the evaluation of the density matrix. Thus if the Green's function
!   $G(ij{\bf k},\omega_n)$ has been determined numerically over all Fermionic
!   Matsubara frequencies up to $\pm \omega_{N_{\rm F}}$, then the density
!   matrix is approximated by
!   $$ \gamma_{ij{\bf k}}=\frac{1}{\beta}\sum_{n\;{\rm odd}}^{\pm N_{\rm F}}
!    \left[G(ij{\bf k},\omega_n)+\frac{a_2}{\omega_n^2}
!    -\frac{a_4}{\omega_n^4}\right]+\frac{1}{\beta}\left[\frac{\beta}{2}a_1
!    -\frac{\beta^2}{4}a_2+\frac{\beta^4}{48}a_4\right], $$
!   where $a_1$, $a_2$ and $a_4$ are chosen so that the Green's function is
!   equal to
!   $$ g(z)=\frac{a_1}{z}+\frac{a_2}{z^2}+\frac{a_3}{z^3}+\frac{a_4}{z^4} $$
!   at the points $n\in \{-n_{\rm F}, -n_{\rm F}+2, n_{\rm F}-2, n_{\rm F}\}$.
!
! !REVISION HISTORY:
!   Created April 2018 (A. Davydov)
!   Increased Laurent series order to 4, December 2023 (JKD)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(in) :: ge(4)
! local variables
integer iw
real(8) sm2,sm4,t1,t2,t3
real(8) c2,c3,c4,d2,d3,d4
complex(8) a1,a2,a4,z1,z2
! π / β
t1=pi*kboltz*tempk
t2=1.d0/(8.d0*dble(nwfm-1))
t3=dble(nwfm)
c2=t3**2
c3=c2*t3
c4=c3*t3
t3=dble(nwfm-2)
d2=t3**2
d3=d2*t3
d4=d3*t3
! determine the coefficients a1, a2 and a4
a1=(d3*(ge(2)-ge(3))+c3*(ge(4)-ge(1)))*t1*t2
a1=cmplx(-aimag(a1),dble(a1),8)
z1=ge(2)+ge(3)
z2=ge(1)+ge(4)
t3=t1**2
a2=(d4*z1-c4*z2)*t3*t2
t3=t3**2
a4=(d2*z1-c2*z2)*c2*d2*t3*t2
! evaluate the sums analytically over all Matsubara frequencies
t1=1.d0/(2.d0*kboltz*tempk)
t2=t1**2
gwtails=a1*t1-a2*t2+a4*(t2**2)/3.d0
! remove the contributions over the finite set of frequencies
sm2=0.d0
sm4=0.d0
do iw=0,nwfm
  t1=aimag(wfm(iw))**2
  sm2=sm2+1.d0/t1
  sm4=sm4+1.d0/t1**2
end do
gwtails=gwtails+sm2*a2-sm4*a4
end function
!EOC

