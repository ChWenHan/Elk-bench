
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengvc
use modmain
implicit none
! local variables
integer ig,ifg,n1
integer i1,i2,i3,j1,j2,j3
real(8) gk2
! find optimal grid size for |G| < 2 gkmax
gk2=2.d0*gkmax+epslat
ngdgc(1:3)=int(gk2*sqrt(avec(1,1:3)**2+avec(2,1:3)**2+avec(3,1:3)**2)/pi)+1
! find next largest FFT-compatible grid size
call nfftifc(npfftgc,ngdgc(1))
call nfftifc(npfftgc,ngdgc(2))
call nfftifc(npfftgc,ngdgc(3))
! total number of points in coarse grid
ngtc=ngdgc(1)*ngdgc(2)*ngdgc(3)
! find the number of vectors with |G| < 2 gkmax
ngvc=ngvec
do ig=2,ngvec
  if (gc(ig) > gk2) then
    ngvc=ig-1
    exit
  end if
end do
! number of complex FFT elements for real-complex transforms
n1=ngdgc(1)/2+1
nfgrzc=n1*ngdgc(2)*ngdgc(3)
! Fourier transform index
if (allocated(igfc)) deallocate(igfc)
allocate(igfc(ngvc))
if (allocated(igrzfc)) deallocate(igrzfc)
allocate(igrzfc(nfgrzc))
igrzfc(1:nfgrzc)=ngvc+1
do ig=1,ngvc
  i1=ivg(1,ig)
  i2=ivg(2,ig)
  i3=ivg(3,ig)
  if (i1 >= 0) then
    j1=i1
  else
    j1=ngdgc(1)+i1
  end if
  if (i2 >= 0) then
    j2=i2
  else
    j2=ngdgc(2)+i2
  end if
  if (i3 >= 0) then
    j3=i3
  else
    j3=ngdgc(3)+i3
  end if
  igfc(ig)=j3*ngdgc(2)*ngdgc(1)+j2*ngdgc(1)+j1+1
! map from real-complex FFT index to G-vector index
  if (i1 >= 0) then
    ifg=j3*ngdgc(2)*n1+j2*n1+j1+1
    igrzfc(ifg)=ig
  end if
end do
end subroutine

