
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findsymcrys
! !INTERFACE:
subroutine findsymcrys
! !USES:
use modmain
use modmpi
use modtest
! !DESCRIPTION:
!   Finds the complete set of symmetries which leave the crystal structure
!   (including the magnetic fields) invariant. A crystal symmetry is of the form
!   $\{\alpha_S|\alpha_R|{\bf t}\}$, where ${\bf t}$ is a translation vector,
!   $\alpha_R$ is a spatial rotation operation and $\alpha_S$ is a global spin
!   rotation. Note that the order of operations is important and defined to be
!   from right to left, i.e. translation followed by spatial rotation followed
!   by spin rotation. In the case of spin-orbit coupling $\alpha_S=\alpha_R$. In
!   order to determine the translation vectors, the entire atomic basis is
!   shifted so that the first atom in the smallest set of atoms of the same
!   species is at the origin. Then all displacement vectors between atoms in
!   this set are checked as possible symmetry translations. If the global
!   variable {\tt tshift} is set to {\tt .false.} then the shift is not
!   performed. See L. M. Sandratskii and P. G. Guletskii, {\it J. Phys. F: Met.
!   Phys.} {\bf 16}, L43 (1986) and the routine {\tt findsym}.
!
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!   Modified for trimvg=.false., November 2023 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ia,ja,is,is0
integer isym,nsym,i,n
integer lspl(48),lspn(48),ilspl
real(8) v0(3),v1(3),v2(3),t1
real(8) apl1(3,maxatoms,maxspecies)
real(8) apl2(3,maxatoms,maxspecies)
! automatic arrays
integer iea(natmmax,nspecies,48)
real(8) vtl(3,natmmax**2+1)
! allocate global equivalent atom arrays
if (allocated(ieqatom)) deallocate(ieqatom)
allocate(ieqatom(natmmax,nspecies,maxsymcrys))
if (allocated(eqatoms)) deallocate(eqatoms)
allocate(eqatoms(natmmax,natmmax,nspecies))
! find the smallest set of atoms
is0=1
do is=1,nspecies
  if (natoms(is) < natoms(is0)) is0=is
end do
if (natmtot > 0) then
! position of first atom in the smallest atom set
  v0(:)=atposl(:,1,is0)
! shift basis so that the first atom in the smallest atom set is at the origin
  do is=1,nspecies
    do ia=1,natoms(is)
! shift atom
      apl1(:,ia,is)=atposl(:,ia,is)-v0(:)
! map lattice coordinates back to [0,1)
      call r3frac(epslat,apl1(:,ia,is))
    end do
  end do
else
  v0(:)=0.d0
end if
! determine possible translation vectors from smallest set of atoms
n=1
vtl(:,1)=0.d0
do ia=1,natoms(is0)
  do ja=2,natoms(is0)
! compute difference between two atom vectors
    v1(:)=apl1(:,ia,is0)-apl1(:,ja,is0)
! map lattice coordinates to [0,1)
    call r3frac(epslat,v1)
! check if vector has any component along electric field
    if (tefield) then
      call r3mv(avec,v1,v2)
      t1=efieldc(1)*v2(1)+efieldc(2)*v2(2)+efieldc(3)*v2(3)
      if (abs(t1) > epslat) goto 10
    end if
    do i=1,n
      t1=abs(vtl(1,i)-v1(1))+abs(vtl(2,i)-v1(2))+abs(vtl(3,i)-v1(3))
      if (t1 < epslat) goto 10
    end do
    n=n+1
    vtl(:,n)=v1(:)
10 continue
  end do
end do
! no translations required when symtype=0,2 (F. Cricchio)
if (symtype /= 1) n=1
eqatoms(:,:,:)=.false.
nsymcrys=0
! loop over all possible translations
do i=1,n
! construct new array with translated positions
  do is=1,nspecies
    do ia=1,natoms(is)
      apl2(:,ia,is)=apl1(:,ia,is)+vtl(:,i)
    end do
  end do
! find the symmetries for current translation
  call findsym(apl1,apl2,nsym,lspl,lspn,iea)
  do isym=1,nsym
    nsymcrys=nsymcrys+1
    if (nsymcrys > maxsymcrys) then
      write(*,*)
      write(*,'("Error(findsymcrys): too many crystal symmetries")')
      write(*,'(" Adjust maxsymcrys in modmain and recompile code")')
      write(*,*)
      stop
    end if
    vtlsymc(:,nsymcrys)=vtl(:,i)
    lsplsymc(nsymcrys)=lspl(isym)
    lspnsymc(nsymcrys)=lspn(isym)
    do is=1,nspecies
      do ia=1,natoms(is)
        ja=iea(ia,is,isym)
        ieqatom(ia,is,nsymcrys)=ja
        eqatoms(ia,ja,is)=.true.
        eqatoms(ja,ia,is)=.true.
      end do
    end do
  end do
end do
tsyminv=.false.
do isym=1,nsymcrys
! check if inversion symmetry is present
  i=lsplsymc(isym)
  if (all(symlat(:,:,i) == -symlat(:,:,1))) then
    tsyminv=.true.
! make inversion the second symmetry element (the identity is the first)
    v1(:)=vtlsymc(:,isym); vtlsymc(:,isym)=vtlsymc(:,2); vtlsymc(:,2)=v1(:)
    i=lsplsymc(isym); lsplsymc(isym)=lsplsymc(2); lsplsymc(2)=i
    i=lspnsymc(isym); lspnsymc(isym)=lspnsymc(2); lspnsymc(2)=i
    do is=1,nspecies
      do ia=1,natoms(is)
        i=ieqatom(ia,is,isym)
        ieqatom(ia,is,isym)=ieqatom(ia,is,2)
        ieqatom(ia,is,2)=i
      end do
    end do
    exit
  end if
end do
if (tshift) then
  if (tsyminv) then
! if inversion exists then shift basis so that inversion center is at origin
    v1(:)=v1(:)/2.d0
  else
    v1(:)=0.d0
  end if
else
  v1(:)=v0(:)
end if
do is=1,nspecies
  do ia=1,natoms(is)
! shift atom
    atposl(:,ia,is)=apl1(:,ia,is)+v1(:)
! map lattice coordinates back to [0,1)
    call r3frac(epslat,atposl(:,ia,is))
! map lattice coordinates to [-0.5,0.5) if inversion exists
    if (tsyminv) then
      do i=1,3
        if (atposl(i,ia,is) > 0.5d0) atposl(i,ia,is)=atposl(i,ia,is)-1.d0
      end do
    end if
! determine the new Cartesian coordinates
    call r3mv(avec,atposl(:,ia,is),atposc(:,ia,is))
  end do
end do
do isym=1,nsymcrys
! recalculate crystal symmetry translation vectors
  ilspl=isymlat(lsplsymc(isym))
  v2(:)=symlat(:,1,ilspl)*v1(1) &
       +symlat(:,2,ilspl)*v1(2) &
       +symlat(:,3,ilspl)*v1(3)
  vtlsymc(:,isym)=vtlsymc(:,isym)-v1(:)+v2(:)
  call r3frac(epslat,vtlsymc(:,isym))
! translation vector in Cartesian coordinates
  call r3mv(avec,vtlsymc(:,isym),vtcsymc(:,isym))
! set flag for zero translation vector
  t1=abs(vtlsymc(1,isym))+abs(vtlsymc(2,isym))+abs(vtlsymc(3,isym))
  if (t1 < epslat) then
    tv0symc(isym)=.true.
  else
    tv0symc(isym)=.false.
  end if
end do
! check inversion does not include a translation
if (tsyminv) then
  if (.not.tv0symc(2)) tsyminv=.false.
end if
if (tshift.and.(natmtot > 0)) then
  v1(:)=atposl(:,1,is0)-v0(:)
  call r3frac(epslat,v1)
  t1=abs(v1(1))+abs(v1(2))+abs(v1(3))
  if (mp_mpi.and.(t1 > epslat)) then
    write(*,*)
    write(*,'("Info(findsymcrys): atomic basis shift (lattice) :")')
    write(*,'(3G18.10)') v1(:)
    write(*,'("See GEOMETRY.OUT for new atomic positions")')
  end if
end if
! write number of crystal symmetries to test file
call writetest(705,'number of crystal symmetries',iv=nsymcrys)
end subroutine
!EOC

