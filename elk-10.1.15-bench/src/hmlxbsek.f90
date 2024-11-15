
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlxbsek(ik2)
use modmain
implicit none
! arguments
integer, intent(in) :: ik2
! local variables
integer ik1,ist1,ist2,jst1,jst2
integer i1,i2,j1,j2,a1,a2,b1,b2
integer is,ias,l
real(8) t0
complex(8) z1
! automatic arrays
integer ngp(nspnfv)
! allocatable arrays
integer, allocatable :: igpig(:,:)
complex(4), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:)
complex(4), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(4), allocatable :: crhomt(:,:),crhoir(:)
complex(4), allocatable :: cvclmt(:,:,:),cvclir(:,:)
! external functions
complex(8), external :: zcfinp
! allocate local arrays
allocate(igpig(ngkmax,nspnfv))
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv),wfir1(ngtc,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv),wfir2(ngtc,nspinor,nstsv))
allocate(crhomt(npcmtmax,natmtot),crhoir(ngtc))
allocate(cvclmt(npcmtmax,natmtot,nvcbse),cvclir(ngtc,nvcbse))
! calculate the wavefunctions for all states of k-point ik2
call genwfsvp_sp(.false.,.false.,nstsv,[0],ngdgc,igfc,vkl(:,ik2),ngp,igpig, &
 wfmt2,ngtc,wfir2)
l=0
do i2=1,nvbse
  ist2=istbse(i2,ik2)
  do j2=1,ncbse
    jst2=jstbse(j2,ik2)
    a2=ijkbse(i2,j2,ik2)
    l=l+1
! calculate the complex overlap density
    call gencrho(.true.,.true.,ngtc,wfmt2(:,:,:,ist2),wfir2(:,:,ist2), &
     wfmt2(:,:,:,jst2),wfir2(:,:,jst2),crhomt,crhoir)
! compute the Coulomb potential
    call gencvclmt(nrcmt,nrcmti,nrcmtmax,rlcmt,wprcmt,npcmtmax,crhomt, &
     cvclmt(:,:,l))
    call cpotcoul(nrcmt,nrcmti,npcmt,nrcmtmax,rlcmt,ngdgc,igfc,ngvc,gc,gclg, &
     ngvec,jlgrmt,ylmg,sfacg,crhoir,npcmtmax,cvclmt(:,:,l),cvclir(:,l))
    cvclir(:,l)=cvclir(:,l)*cfrc(:)
  end do
end do
t0=occmax*wkptnr
! start loop over ik1
do ik1=1,nkptnr
  if (ik1 == ik2) then
    wfmt1(:,:,:,:)=wfmt2(:,:,:,:)
    wfir1(:,:,:)=wfir2(:,:,:)
  else
    call genwfsvp_sp(.false.,.false.,nstsv,[0],ngdgc,igfc,vkl(:,ik1),ngp,igpig,&
     wfmt1,ngtc,wfir1)
  end if
  do i1=1,nvbse
    ist1=istbse(i1,ik1)
    do j1=1,ncbse
      jst1=jstbse(j1,ik1)
      a1=ijkbse(i1,j1,ik1)
! calculate the complex overlap density
      call gencrho(.true.,.true.,ngtc,wfmt1(:,:,:,ist1),wfir1(:,:,ist1), &
       wfmt1(:,:,:,jst1),wfir1(:,:,jst1),crhomt,crhoir)
      l=0
      do i2=1,nvbse
        ist2=istbse(i2,ik2)
        do j2=1,ncbse
          jst2=jstbse(j2,ik2)
          a2=ijkbse(i2,j2,ik2)
          l=l+1
! compute the matrix element
          z1=t0*zcfinp(crhomt,crhoir,cvclmt(:,:,l),cvclir(:,l))
          hmlbse(a1,a2)=hmlbse(a1,a2)+z1
! compute off-diagonal blocks if required
          if (bsefull) then
            b1=a1+nbbse
            b2=a2+nbbse
            hmlbse(b1,b2)=hmlbse(b1,b2)-conjg(z1)
! conjugate the potential
            do ias=1,natmtot
              is=idxis(ias)
              call cfmtconj(nrcmt(is),nrcmti(is),npcmt(is),cvclmt(:,ias,l))
            end do
            cvclir(:,l)=conjg(cvclir(:,l))
            z1=t0*zcfinp(crhomt,crhoir,cvclmt(:,:,l),cvclir(:,l))
            hmlbse(a1,b2)=hmlbse(a1,b2)+z1
            hmlbse(b1,a2)=hmlbse(b1,a2)-conjg(z1)
          end if
        end do
      end do
    end do
  end do
end do
deallocate(igpig,wfmt1,wfmt2,wfir1,wfir2)
deallocate(crhomt,crhoir,cvclmt,cvclir)
end subroutine

