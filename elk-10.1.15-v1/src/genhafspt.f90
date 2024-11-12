
! Copyright (C) 2023 E. Harris-Lee, S. Sharma and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genhafspt(evecsv,pmat,h)
use modmain
use modtddft
implicit none
! arguments
complex(8), intent(in) :: evecsv(nstsv,nstsv),pmat(nstsv,nstsv,3)
complex(8), intent(inout) :: h(nstsv,nstsv)
! local variables
integer i,j,l
real(8) ca,v(3)
complex(8) zm(2,2)
! allocatable arrays
complex(8), allocatable :: smat(:,:,:,:)
allocate(smat(nstsv,nstsv,2,2))
! generate the spin operator matrix elements in the second-variational basis
call gensmatk(evecsv,smat)
! coupling constant of the external spin-polarised A-field (-1/c)
ca=-1.d0/solsc
! add the spin-current operator to the Hamiltonian
do l=1,3
  v(:)=ca*afspt(l,:,itimes)
  zm(1,1)=v(3)
  zm(1,2)=cmplx(v(1),-v(2),8)
  zm(2,1)=cmplx(v(1),v(2),8)
  zm(2,2)=-v(3)
  do j=1,2
    do i=1,2
      call zhemm('L','U',nstsv,nstsv,zm(i,j),pmat(:,:,l),nstsv,smat(:,:,i,j), &
       nstsv,zone,h,nstsv)
    end do
  end do
end do
deallocate(smat)
end subroutine

