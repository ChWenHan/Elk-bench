
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gencrho(tsh,tspc,ngt,wfmt1,wfir1,wfmt2,wfir2,crhomt,crhoir)
use modmain
implicit none
! arguments
logical, intent(in) :: tsh,tspc
integer, intent(in) :: ngt
complex(4), intent(in) :: wfmt1(npcmtmax,natmtot,*),wfir1(ngt,*)
complex(4), intent(in) :: wfmt2(npcmtmax,natmtot,*),wfir2(ngt,*)
complex(4), intent(out) :: crhomt(npcmtmax,natmtot),crhoir(ngt)
! local variables
integer is,ias
! muffin-tin part
do ias=1,natmtot
  is=idxis(ias)
  if (tsh) then
    if (tspc.and.spinpol) then
! contract over spin
      call crho2(npcmt(is),wfmt1(:,ias,1),wfmt1(:,ias,2),wfmt2(:,ias,1), &
       wfmt2(:,ias,2),crhomt(:,ias))
    else
! no spin contraction
      call crho1(npcmt(is),wfmt1(:,ias,1),wfmt2(:,ias,1),crhomt(:,ias))
    end if
! convert to spherical harmonics
    call cfshtip(nrcmt(is),nrcmti(is),crhomt(:,ias))
  else
    if (tspc.and.spinpol) then
      call crho2(npcmt(is),wfmt1(:,ias,1),wfmt1(:,ias,2),wfmt2(:,ias,1), &
       wfmt2(:,ias,2),crhomt(:,ias))
    else
      call crho1(npcmt(is),wfmt1(:,ias,1),wfmt2(:,ias,1),crhomt(:,ias))
    end if
  end if
end do
! interstitial part
if (tspc.and.spinpol) then
  call crho2(ngt,wfir1,wfir1(:,2),wfir2,wfir2(:,2),crhoir)
else
  call crho1(ngt,wfir1,wfir2,crhoir)
end if
return

contains

pure subroutine crho1(n,wf1,wf2,crho)
implicit none
integer, intent(in) :: n
complex(4), intent(in) :: wf1(n),wf2(n)
complex(4), intent(out) :: crho(n)
crho(1:n)=conjg(wf1(1:n))*wf2(1:n)
end subroutine

pure subroutine crho2(n,wf11,wf12,wf21,wf22,crho)
implicit none
integer, intent(in) :: n
complex(4), intent(in) :: wf11(n),wf12(n),wf21(n),wf22(n)
complex(4), intent(out) :: crho(n)
crho(1:n)=conjg(wf11(1:n))*wf21(1:n)+conjg(wf12(1:n))*wf22(1:n)
end subroutine

end subroutine

