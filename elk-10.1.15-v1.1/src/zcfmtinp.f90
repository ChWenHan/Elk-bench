
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure complex(8) function zcfmtinp(nr,nri,wr,cfmt1,cfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr)
complex(4), intent(in) :: cfmt1(*),cfmt2(*)
! local variables
integer n,ir,i
! compute the dot-products for each radial point and integrate over r
zcfmtinp=0.d0
if (lmaxi == 1) then
  do ir=1,nri
    i=4*(ir-1)+1
    zcfmtinp=zcfmtinp+wr(ir) &
     *(conjg(cfmt1(i))*cfmt2(i) &
      +conjg(cfmt1(i+1))*cfmt2(i+1) &
      +conjg(cfmt1(i+2))*cfmt2(i+2) &
      +conjg(cfmt1(i+3))*cfmt2(i+3))
  end do
  i=4*nri+1
else
  i=1
  n=lmmaxi-1
  do ir=1,nri
    zcfmtinp=zcfmtinp+wr(ir)*dot_product(cfmt1(i:i+n),cfmt2(i:i+n))
    i=i+lmmaxi
  end do
end if
n=lmmaxo-1
do ir=nri+1,nr
  zcfmtinp=zcfmtinp+wr(ir)*dot_product(cfmt1(i:i+n),cfmt2(i:i+n))
  i=i+lmmaxo
end do
end function

