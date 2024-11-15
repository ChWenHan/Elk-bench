
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmatk(tspndg,tlmdg,lmin,lmax,ias,nst,idx,ngp,apwalm,evecfv, &
 evecsv,ld,dmat)
use modmain
implicit none
! arguments
logical, intent(in) :: tspndg,tlmdg
integer, intent(in) :: lmin,lmax,ias
integer, intent(in) :: nst,idx(*),ngp(nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
integer, intent(in) :: ld
complex(8), intent(out) :: dmat(ld,nspinor,ld,nspinor,nst)
! local variables
integer ispn,jspn,ist,is
integer nrc,nrci,irco
integer l,lma,lmb,lm1,lm2
integer npci,ni,no,i1,i2
complex(8) zsm
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:)
if (lmin < 0) then
  write(*,*)
  write(*,'("Error(gendmatk): lmin < 0 : ",I8)') lmin
  write(*,*)
  stop
end if
if (lmax > lmaxo) then
  write(*,*)
  write(*,'("Error(gendmatk): lmax > lmaxo : ",2I8)') lmax,lmaxo
  write(*,*)
  stop
end if
is=idxis(ias)
nrc=nrcmt(is)
nrci=nrcmti(is)
irco=nrci+1
npci=npcmti(is)
ni=npci-1
no=npcmt(is)-npci-1
! generate the second-variational wavefunctions
allocate(wfmt(npcmtmax,nspinor,nst))
call wfmtsv(.true.,lradstp,is,ias,nst,idx,ngp,apwalm,evecfv,evecsv,npcmtmax, &
 wfmt)
! zero the density matrix
dmat(:,:,:,:,:)=0.d0
! loop over second-variational states
do ist=1,nst
  do ispn=1,nspinor
    do jspn=1,nspinor
      if (tspndg.and.(ispn /= jspn)) cycle
      do l=lmin,lmax
        lma=l**2+1; lmb=lma+2*l
        do lm1=lma,lmb
          do lm2=lma,lmb
            if (tlmdg.and.(lm1 /= lm2)) cycle
            if (l <= lmaxi) then
              zsm=sum(wfmt(lm1:lm1+ni:lmmaxi,ispn,ist) &
               *conjg(wfmt(lm2:lm2+ni:lmmaxi,jspn,ist))*wrcmt(1:nrci,is))
            else
              zsm=0.d0
            end if
            i1=npci+lm1; i2=npci+lm2
            zsm=zsm+sum(wfmt(i1:i1+no:lmmaxo,ispn,ist) &
             *conjg(wfmt(i2:i2+no:lmmaxo,jspn,ist))*wrcmt(irco:nrc,is))
            dmat(lm1,ispn,lm2,jspn,ist)=zsm
          end do
        end do
      end do
    end do
  end do
! end loop over second-variational states
end do
deallocate(wfmt)
end subroutine

