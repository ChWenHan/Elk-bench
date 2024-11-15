
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynrtoq(vpl,dr,dq)
use modmain
use modphonon
implicit none
! arguments
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: dr(nbph,nbph,nqptnr)
complex(8), intent(out) :: dq(nbph,nbph)
! local variables
integer isym,lspl,ilspl
integer i1,i2,i3,ir
real(8) s(3,3),v(3),t1
complex(8) z1
! automatic arrays
complex(8) d(nbph,nbph)
dq(:,:)=0.d0
! loop over crystal symmetries
do isym=1,nsymcrys
! index to spatial rotation in lattice point group
  lspl=lsplsymc(isym)
! the inverse of the spatial symmetry
  ilspl=isymlat(lspl)
! symmetry matrix in lattice coordinates
  s(:,:)=dble(symlat(:,:,ilspl))
! operate with inverse symmetry matrix on vpl
  call r3mtv(s,vpl,v)
! construct dynamical matrix for rotated vpl
  d(:,:)=0.d0
! loop over R-vectors
  ir=0
  do i3=ngridq(3)/2-ngridq(3)+1,ngridq(3)/2
    do i2=ngridq(2)/2-ngridq(2)+1,ngridq(2)/2
      do i1=ngridq(1)/2-ngridq(1)+1,ngridq(1)/2
        ir=ir+1
        t1=-twopi*(v(1)*dble(i1)+v(2)*dble(i2)+v(3)*dble(i3))
        z1=cmplx(cos(t1),sin(t1),8)
        d(:,:)=d(:,:)+z1*dr(:,:,ir)
      end do
    end do
  end do
! apply symmetry operation to dynamical matrix and add to total
  call dynsymapp(isym,v,d,dq)
! end loop over symmetries
end do
! normalise by the number of symmetry operations
t1=1.d0/dble(nsymcrys)
dq(:,:)=t1*dq(:,:)
! add the non-analytic term if required
if (tphnat) call dynqnat(1,vpl,dq)
end subroutine

