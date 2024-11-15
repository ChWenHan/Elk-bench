
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genffacgp
! !INTERFACE:
pure subroutine genffacgp(ngp,gpc,ld,ffacgp)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   gpc    : length of G+p-vectors (in,real(ngp))
!   ld     : leading dimension (in,integer)
!   ffacgp : form factors (out,real(ld,nspecies))
! !DESCRIPTION:
!   Generates the form factors used to determine the smooth characteristic
!   function. See {\tt gencfun} for details.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: gpc(ngp)
integer, intent(in) :: ld
real(8), intent(out) :: ffacgp(ld,nspecies)
! local variables
integer is,ig
real(8) r,t0,t1
t0=fourpi/omega
do is=1,nspecies
  r=rmt(is)
  do ig=1,ngp
    if (gpc(ig) > epslat) then
      t1=gpc(ig)*r
      ffacgp(ig,is)=t0*(sin(t1)-t1*cos(t1))/(gpc(ig)**3)
    else
      ffacgp(ig,is)=(t0/3.d0)*(r**3)
    end if
  end do
end do
end subroutine
!EOC

