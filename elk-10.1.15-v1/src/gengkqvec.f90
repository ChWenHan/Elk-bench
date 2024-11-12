
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengkqvec(iq)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq
! local variables
integer ik,jspn
real(8) vl(3),vc(3)
! loop over non-reduced k-point set
do ik=1,nkptnr
! k+q-vectors in lattice and Cartesian coordinates
  vkql(1:3,ik)=vkl(1:3,ik)+vql(1:3,iq)
  vkqc(1:3,ik)=vkc(1:3,ik)+vqc(1:3,iq)
  do jspn=1,nspnfv
    vl(1:3)=vkql(1:3,ik)
    vc(1:3)=vkqc(1:3,ik)
! spin-spiral case
    if (spinsprl) then
      if (jspn == 1) then
        vl(1:3)=vl(1:3)+0.5d0*vqlss(1:3)
        vc(1:3)=vc(1:3)+0.5d0*vqcss(1:3)
      else
        vl(1:3)=vl(1:3)-0.5d0*vqlss(1:3)
        vc(1:3)=vc(1:3)-0.5d0*vqcss(1:3)
      end if
    end if
! generate the G+k+q-vectors
    call gengkvec(ngvc,ivg,vgc,vl,vc,gkmax,ngkmax,ngkq(jspn,ik), &
     igkqig(:,jspn,ik),vgkql(:,:,jspn,ik),vgkqc(:,:,jspn,ik),gkqc(:,jspn,ik))
! generate structure factors for the G+k+q-vectors
    call gensfacgp(ngkq(jspn,ik),vgkqc(:,:,jspn,ik),ngkmax,sfacgkq(:,:,jspn,ik))
  end do
end do
end subroutine

