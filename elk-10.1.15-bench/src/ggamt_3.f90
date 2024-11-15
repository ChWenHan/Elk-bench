
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ggamt_3(is,np,vx,vc,wx,wc,dtdg2r)
use modmain
implicit none
! arguments
integer, intent(in) :: is,np
real(8), intent(inout) :: vx(np),vc(np)
real(8), intent(in) :: wx(np),wc(np)
real(8), intent(in) :: dtdg2r(np)
! local variables
integer nr,nri
! automatic arrays
real(8) rfmt1(np),rfmt2(np)
nr=nrmt(is)
nri=nrmti(is)
!------------------!
!     exchange     !
!------------------!
rfmt1(1:np)=wx(1:np)*dtdg2r(1:np)
call rfsht(nr,nri,rfmt1,rfmt2)
call grad2rfmt(nr,nri,rlmt(:,-1,is),rlmt(:,-2,is),wcrmt(:,:,is),rfmt2,rfmt1)
call rbsht(nr,nri,rfmt1,rfmt2)
vx(1:np)=vx(1:np)+rfmt2(1:np)
!---------------------!
!     correlation     !
!---------------------!
rfmt1(1:np)=wc(1:np)*dtdg2r(1:np)
call rfsht(nr,nri,rfmt1,rfmt2)
call grad2rfmt(nr,nri,rlmt(:,-1,is),rlmt(:,-2,is),wcrmt(:,:,is),rfmt2,rfmt1)
call rbsht(nr,nri,rfmt1,rfmt2)
vc(1:np)=vc(1:np)+rfmt2(1:np)
end subroutine

