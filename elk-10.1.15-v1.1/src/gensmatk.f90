
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gensmatk(evecsv,smat)
use modmain
implicit none
! arguments
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: smat(nstsv,nstsv,2,2)
! local variables
integer ist,jst,i
! external functions
complex(8), external :: zdotc
i=nstfv+1
do jst=1,nstsv
  do ist=1,jst
    smat(ist,jst,1,1)=zdotc(nstfv,evecsv(1,ist),1,evecsv(1,jst),1)
    smat(ist,jst,2,1)=zdotc(nstfv,evecsv(i,ist),1,evecsv(1,jst),1)
    smat(ist,jst,1,2)=zdotc(nstfv,evecsv(1,ist),1,evecsv(i,jst),1)
    smat(ist,jst,2,2)=zdotc(nstfv,evecsv(i,ist),1,evecsv(i,jst),1)
  end do
end do
! set the lower triangular parts
do jst=1,nstsv
  do ist=1,jst-1
    smat(jst,ist,1,1)=conjg(smat(ist,jst,1,1))
    smat(jst,ist,2,1)=conjg(smat(ist,jst,1,2))
    smat(jst,ist,1,2)=conjg(smat(ist,jst,2,1))
    smat(jst,ist,2,2)=conjg(smat(ist,jst,2,2))
  end do
end do
end subroutine

