
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine straingkq
use modmain
use modulr
implicit none
integer is,ia,ig
integer nppt,ik,igk
integer jspn,iq
real(8) ta(3,3),tb(3,3),vc(3)
if ((istrain < 1).or.(istrain > nstrain)) return
! compute the strained lattice vectors
avec(1:3,1:3)=avec0(1:3,1:3)+deltast*strain(1:3,1:3,istrain)
! generate the strained reciprocal lattice vectors and unit cell volume
call reciplat(avec,bvec,omega,omegabz)
! determine the transformation matrix to the strained vectors
call r3mm(avec,ainv,ta)
call r3mm(bvec,binv,tb)
! recalculate all required variables which depend on avec
call r3minv(avec,ainv)
call r3minv(bvec,binv)
call r3mv(bvec,vqlss,vqcss)
do is=1,nspecies
  do ia=1,natoms(is)
    call r3mv(avec,atposl(:,ia,is),atposc(:,ia,is))
  end do
end do
call r3mv(bvec,vecql,vecqc)
call r3mv(ainv,efieldc,efieldl)
call r3mv(ainv,afieldc,afieldl)
call symmetry
! apply the transformation matrix to the G-vectors
do ig=1,ngtot
  vc(1:3)=vgc(1:3,ig)
  call r3mv(tb,vc,vgc(:,ig))
  gc(ig)=sqrt(vgc(1,ig)**2+vgc(2,ig)**2+vgc(3,ig)**2)
end do
! recalculate variables which depend on the G-vectors
call gengclg
call genjlgprmt(lnpsd,ngvec,gc,ngvec,jlgrmt)
call genylmg
call gensfacgp(ngvec,vgc,ngvec,sfacg)
call genffacgp(ngtot,gc,ngtot,ffacg)
call gencfun
call gencfrc
call energynn
! apply the transformation to the k-vectors
do ik=1,nkptnr
  vc(1:3)=vkc(1:3,ik)
  call r3mv(tb,vc,vkc(:,ik))
end do
! apply the transformation to G+k-vectors and recalculate dependent variables
if (xctype(1) < 0) then
  nppt=nkptnr
else
  nppt=nkpt
end if
do ik=1,nppt
  do jspn=1,nspnfv
    do igk=1,ngk(jspn,ik)
      vc(1:3)=vgkc(1:3,igk,jspn,ik)
      call r3mv(tb,vc,vgkc(:,igk,jspn,ik))
      gkc(igk,jspn,ik)=sqrt(vgkc(1,igk,jspn,ik)**2 &
                           +vgkc(2,igk,jspn,ik)**2 &
                           +vgkc(3,igk,jspn,ik)**2)
    end do
    call gensfacgp(ngk(jspn,ik),vgkc(:,:,jspn,ik),ngkmax,sfacgk(:,:,jspn,ik))
  end do
end do
! apply the transformation to the q-vectors if required
if (xctype(1) < 0) then
  do iq=1,nqptnr
    vc(1:3)=vqc(1:3,iq)
    call r3mv(tb,vc,vqc(:,iq))
  end do
  call gengclq
end if
end subroutine

