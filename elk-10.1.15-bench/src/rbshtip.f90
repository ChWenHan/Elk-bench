
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rbshtip(nr,nri,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(inout) :: rfmt(*)
! local variables
integer nro,npi,npo,i
! automatic arrays
real(8) f(max(lmmaxi*nri,lmmaxo*(nr-nri)))
! transform the inner part of the muffin-tin function in-place
npi=lmmaxi*nri
f(1:npi)=rfmt(1:npi)
call dgemm('N','N',lmmaxi,nri,lmmaxi,1.d0,rbshti,lmmaxi,f,lmmaxi,0.d0,rfmt, &
 lmmaxi)
! transform the outer part of the muffin-tin function in-place
nro=nr-nri
npo=lmmaxo*nro
i=npi+1
f(1:npo)=rfmt(i:npi+npo)
call dgemm('N','N',lmmaxo,nro,lmmaxo,1.d0,rbshto,lmmaxo,f,lmmaxo,0.d0,rfmt(i), &
 lmmaxo)
end subroutine

