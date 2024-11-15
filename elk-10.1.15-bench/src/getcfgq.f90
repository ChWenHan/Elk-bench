
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getcfgq(fname,vpl,ng,m,cf)
use modmain
implicit none
! arguments
character(*), intent(in) :: fname
real(8), intent(in) :: vpl(3)
integer, intent(in) :: ng,m
complex(8), intent(out) :: cf(ng,ng,m)
! local variables
integer isym,iq,lspl,ilspl
integer igq,jgq,igp,jgp,i
integer recl,ng_,m_
real(8) vql_(3),si(3,3)
real(8) vgql(3),v(3),t1
complex(8) z1
! automatic arrays
logical done(ng)
integer map(ng)
! allocatable arrays
real(8), allocatable :: vgpl(:,:)
complex(8), allocatable :: cf_(:,:,:),x(:)
! find the equivalent reduced q-point and symmetry which rotates vql to vpl
call findqpt(vpl,isym,iq)
! find the record length
inquire(iolength=recl) vql_,ng_,m_,cf
!$OMP CRITICAL(u245)
open(245,file=trim(fname),form='UNFORMATTED',access='DIRECT',recl=recl)
read(245,rec=iq) vql_,ng_,m_,cf
close(245)
!$OMP END CRITICAL(u245)
t1=abs(vql(1,iq)-vql_(1))+abs(vql(2,iq)-vql_(2))+abs(vql(3,iq)-vql_(3))
if (t1 > epslat) then
  write(*,*)
  write(*,'("Error(getcfgq): differing vectors for q-point ",I8)') iq
  write(*,'(" current : ",3G18.10)') vql(:,iq)
  write(*,'(" file    : ",3G18.10)') vql_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
if (ng /= ng_) then
  write(*,*)
  write(*,'("Error(getcfgq): differing ng for q-point ",I8)') iq
  write(*,'(" current : ",I8)') ng
  write(*,'(" file    : ",I8)') ng_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
if (m /= m_) then
  write(*,*)
  write(*,'("Error(getcfgq): differing m for q-point ",I8)') iq
  write(*,'(" current : ",I8)') m
  write(*,'(" file    : ",I8)') m_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
! if p = q then return
t1=abs(vpl(1)-vql(1,iq))+abs(vpl(2)-vql(2,iq))+abs(vpl(3)-vql(3,iq))
if (t1 < epslat) return
! allocate local arrays
allocate(vgpl(3,ng),cf_(ng,ng,m),x(ng))
! perform translation operation and store in temporary array
if (tv0symc(isym)) then
! translation vector is zero
  cf_(1:ng,1:ng,1:m)=cf(1:ng,1:ng,1:m)
else
! non-zero translation vector gives a phase factor
  v(1:3)=vtcsymc(1:3,isym)
  do igq=1,ng
    t1=-(vgc(1,igq)*v(1)+vgc(2,igq)*v(2)+vgc(3,igq)*v(3))
    x(igq)=cmplx(cos(t1),sin(t1),8)
  end do
  do jgq=1,ng
    z1=conjg(x(jgq))
    do igq=1,ng
      cf_(igq,jgq,1:m)=z1*x(igq)*cf(igq,jgq,1:m)
    end do
  end do
end if
! index to spatial rotation in lattice point group
lspl=lsplsymc(isym)
! the inverse of the spatial symmetry
ilspl=isymlat(lspl)
si(1:3,1:3)=dble(symlat(1:3,1:3,ilspl))
! find the map from {G+q} to {G+p}
map(1:ng)=0
do igp=1,ng
  vgpl(1:3,igp)=dble(ivg(1:3,igp))+vpl(1:3)
end do
done(1:ng)=.false.
i=1
do igq=1,ng
  vgql(1:3)=dble(ivg(1:3,igq))+vql(1:3,iq)
  call r3mtv(si,vgql,v)
  do igp=i,ng
    if (done(igp)) cycle
    t1=abs(v(1)-vgpl(1,igp))+abs(v(2)-vgpl(2,igp))+abs(v(3)-vgpl(3,igp))
    if (t1 < epslat) then
      map(igp)=igq
      done(igp)=.true.
      exit
    end if
  end do
  do igp=i,ng
    if (.not.done(igp)) then
      i=igp
      exit
    end if
  end do
end do
! rotate correlation function (passive transformation)
do jgp=1,ng
  jgq=map(jgp)
  do igp=1,ng
    igq=map(igp)
    if ((igq == 0).or.(jgq == 0)) then
      cf(igp,jgp,1:m)=0.d0
    else
      cf(igp,jgp,1:m)=cf_(igq,jgq,1:m)
    end if
  end do
end do
deallocate(vgpl,cf_,x)
end subroutine

