
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vmatmtsc
use modmain
use moddftu
implicit none
! local variables
integer ispn,jspn,ias,lm
! automatic arrays
complex(8) a(lmmaxo,lmmaxo),b(lmmaxo,lmmaxo),c(lmmaxo,lmmaxo)
lm=min(lmmaxi,lmmaxdm)
! transform the muffin-tin potential matrix elements from a spherical harmonics
! basis to a mixed spherical harmonics/coordinates basis and store in global
! single-precision arrays
do ias=1,natmtot
  if (any(tvmmt(0:lmaxdm,ias))) then
    do jspn=1,nspinor
      do ispn=1,nspinor
! inner part of muffin-tin
        a(1:lmmaxi,1:lmmaxi)=0.d0
        a(1:lm,1:lm)=vmatmt(1:lm,ispn,1:lm,jspn,ias)
        call zgemm('N','N',lmmaxi,lmmaxi,lmmaxi,zone,a,lmmaxo,zfshti,lmmaxi, &
         zzero,b,lmmaxo)
        call zgemm('N','N',lmmaxi,lmmaxi,lmmaxi,zone,zbshti,lmmaxi,b,lmmaxo, &
         zzero,c,lmmaxo)
        vmatmti(1:lmmaxi,1:lmmaxi,ispn,jspn,ias)=c(1:lmmaxi,1:lmmaxi)
! outer part of muffin-tin
        a(:,:)=0.d0
        a(1:lmmaxdm,1:lmmaxdm)=vmatmt(1:lmmaxdm,ispn,1:lmmaxdm,jspn,ias)
        call zgemm('N','N',lmmaxo,lmmaxo,lmmaxo,zone,a,lmmaxo,zfshto,lmmaxo, &
         zzero,b,lmmaxo)
        call zgemm('N','N',lmmaxo,lmmaxo,lmmaxo,zone,zbshto,lmmaxo,b,lmmaxo, &
         zzero,c,lmmaxo)
        vmatmto(:,:,ispn,jspn,ias)=c(:,:)
      end do
    end do
  end if
end do
end subroutine

