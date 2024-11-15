
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvfxcg(gclgq,nm,vfxc)
use modmain
implicit none
! arguments
real(8), intent(in) :: gclgq(ngrf)
integer, intent(in) :: nm
complex(8), intent(out) :: vfxc(nm,nm,nwrf)
! local variables
integer ig,jg,kg,i1,i2,i3
complex(8) z1
! allocatable arrays
real(8), allocatable :: fxcmt(:,:),fxcir(:)
complex(8), allocatable :: fxcg(:)
allocate(fxcmt(npmtmax,natmtot),fxcir(ngtot))
allocate(fxcg(ngvec))
! generate the kernel f_xc in real-space
call genfxcr(.true.,fxcmt,fxcir)
! Fourier transform the kernel to G-space
call zftrf(ngvec,ivg,vgc,fxcmt,fxcir,fxcg)
do ig=1,ngrf
  do jg=1,ngrf
    i1=ivg(1,ig)-ivg(1,jg)
    i2=ivg(2,ig)-ivg(2,jg)
    i3=ivg(3,ig)-ivg(3,jg)
    if ((i1 >= intgv(1,1)).and.(i1 <= intgv(2,1)).and. &
        (i2 >= intgv(1,2)).and.(i2 <= intgv(2,2)).and. &
        (i3 >= intgv(1,3)).and.(i3 <= intgv(2,3))) then
      kg=ivgig(i1,i2,i3)
      if (kg <= ngvec) then
        z1=fxcg(kg)/(gclgq(ig)*gclgq(jg))
        vfxc(ig,jg,1:nwrf)=z1
      else
        vfxc(ig,jg,1:nwrf)=0.d0
      end if
    end if
  end do
end do
deallocate(fxcmt,fxcir,fxcg)
end subroutine

