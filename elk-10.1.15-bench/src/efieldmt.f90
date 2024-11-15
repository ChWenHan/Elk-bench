
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: efieldmt
! !INTERFACE:
subroutine efieldmt
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the average electric field in each muffin-tin from the gradient
!   of the Coulomb potential:
!   \begin{align*}
!    {\bf E}_{\alpha}&\equiv-\frac{3}{4\pi R_{\alpha}^3}
!    \int_{{\rm MT}_\alpha} \nabla V_{\rm C}({\bf r})\,d^3r \\
!    &=-\frac{3}{4\pi R_{\alpha}^3}\int_{{\rm MT}_\alpha}
!    V_{\rm C}({\bf r})\,\hat{\bf n}\,dS,
!   \end{align*}
!   where $R_{\alpha}$ is the radius of muffin-tin $\alpha$.
!
! !REVISION HISTORY:
!   Created April 2024 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ias,nr,nri,i
real(8) t1,t2
! automatic arrays
real(8) grfmt(npmtmax,3)
! external functions
real(8), external :: rfmtint
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),vclmt(:,ias),npmtmax,grfmt)
  t1=-(fourpi/3.d0)*rmt(is)**3
  do i=1,3
    t2=rfmtint(nr,nri,wrmt(:,is),grfmt(:,i))
    efcmt(i,ias)=t2/t1
  end do
end do
end subroutine
!EOC

