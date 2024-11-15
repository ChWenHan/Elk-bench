
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfshtip(nr,nri,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(inout) :: zfmt(*)
! local variables
integer nro,npi,npo,i
! automatic arrays
complex(8) f(max(lmmaxi*nri,lmmaxo*(nr-nri)))
! transform the inner part of the muffin-tin function in-place
npi=lmmaxi*nri
f(1:npi)=zfmt(1:npi)
call zgemm('N','N',lmmaxi,nri,lmmaxi,zone,zfshti,lmmaxi,f,lmmaxi,zzero,zfmt, &
 lmmaxi)
! transform the outer part of the muffin-tin function in-place
nro=nr-nri
npo=lmmaxo*nro
i=npi+1
f(1:npo)=zfmt(i:npi+npo)
call zgemm('N','N',lmmaxo,nro,lmmaxo,zone,zfshto,lmmaxo,f,lmmaxo,zzero,zfmt(i),&
 lmmaxo)
end subroutine

