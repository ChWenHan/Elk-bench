
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine cfshtip(nr,nri,cfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(4), intent(inout) :: cfmt(*)
! local variables
integer nro,npi,npo,i
! automatic arrays
complex(4) f(max(lmmaxi*nri,lmmaxo*(nr-nri)))
! transform the inner part of the muffin-tin function in-place
npi=lmmaxi*nri
f(1:npi)=cfmt(1:npi)
call cgemm('N','N',lmmaxi,nri,lmmaxi,cone,cfshti,lmmaxi,f,lmmaxi,czero,cfmt, &
 lmmaxi)
! transform the outer part of the muffin-tin function in-place
nro=nr-nri
npo=lmmaxo*nro
i=npi+1
f(1:npo)=cfmt(i:npi+npo)
call cgemm('N','N',lmmaxo,nro,lmmaxo,cone,cfshto,lmmaxo,f,lmmaxo,czero,cfmt(i),&
 lmmaxo)
end subroutine

