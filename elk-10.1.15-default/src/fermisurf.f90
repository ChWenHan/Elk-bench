
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine fermisurf
use modmain
use modomp
implicit none
! local variables
integer ik,ist,np,nfp,i
integer ist0,ist1,jst0,jst1
integer l,i1,i2,i3,nthd
real(8) e0,e1,v(3),x,fn
! allocatable arrays
real(8), allocatable :: evalfv(:,:),vpc(:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
! external functions
real(8), external :: sdelta
! initialise universal variables
call init0
! no k-point reduction for the collinear case
reducek0=reducek
if (ndmag == 1) reducek=0
call init1
! read density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! find the new linearisation energies
call linengy
! generate the APW and local-orbital radial functions and integrals
call genapwlofr
! generate the spin-orbit coupling radial functions
call gensocfr
! begin parallel loop over reduced k-points set
call holdthd(nkpt,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv) &
!$OMP NUM_THREADS(nthd)
allocate(evalfv(nstfv,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
!$OMP DO
do ik=1,nkpt
!$OMP CRITICAL(fermisurf_)
  write(*,'("Info(fermisurf): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(fermisurf_)
! solve the first- and second-variational eigenvalue equations
  call eveqn(ik,evalfv,evecfv,evecsv)
! end loop over reduced k-points set
end do
!$OMP END DO
deallocate(evalfv,evecfv,evecsv)
!$OMP END PARALLEL
call freethd(nthd)
! if iterative diagonalisation is used the eigenvalues must be reordered
if (tefvit.and.(.not.spinpol)) then
  do ik=1,nkpt
    call sort(nstsv,evalsv(:,ik))
  end do
end if
! generate the plotting point grid in Cartesian coordinates (this has the same
! arrangement as the k-point grid)
np=np3d(1)*np3d(2)*np3d(3)
allocate(vpc(3,np))
call plotpt3d(vpc)
do i=1,np
  v(:)=vpc(:,i)
  call r3mv(bvec,v,vpc(:,i))
end do
! number of files to plot (2 for collinear magnetism, 1 otherwise)
if (ndmag == 1) then
  nfp=2
else
  nfp=1
end if
do l=1,nfp
  if (nfp == 2) then
    if (l == 1) then
      open(50,file='FERMISURF_UP.OUT',form='FORMATTED',action='WRITE')
      jst0=1; jst1=nstfv
    else
      open(50,file='FERMISURF_DN.OUT',form='FORMATTED',action='WRITE')
      jst0=nstfv+1; jst1=2*nstfv
    end if
  else
    open(50,file='FERMISURF.OUT',form='FORMATTED',action='WRITE')
    jst0=1; jst1=nstsv
  end if
  if ((task == 100).or.(task == 101)) then
! find the range of eigenvalues which contribute to the Fermi surface (Lars)
    ist0=jst1; ist1=jst0
    do ist=jst0,jst1
      e0=minval(evalsv(ist,:)); e1=maxval(evalsv(ist,:))
! determine if the band crosses the Fermi energy
      if ((e0 < efermi).and.(e1 > efermi)) then
        ist0=min(ist0,ist); ist1=max(ist1,ist)
      end if
    end do
  else
! use all eigenvalues
    ist0=jst0; ist1=jst1
  end if
  if ((task == 100).or.(task == 103)) then
    write(50,'(3I6," : grid size")') np3d(:)
    i=0
    do i3=0,ngridk(3)-1
      do i2=0,ngridk(2)-1
        do i1=0,ngridk(1)-1
          i=i+1
          ik=ivkik(i1,i2,i3)
          if (task == 100) then
! write product of eigenstates minus the Fermi energy
            fn=product(evalsv(ist0:ist1,ik)-efermi)
          else
! write single smooth delta function at the Fermi energy
            fn=0.d0
            do ist=ist0,ist1
              x=(evalsv(ist,ik)-efermi)/swidth
              fn=fn+sdelta(stype,x)/swidth
            end do
          end if
          write(50,'(4G18.10)') vpc(:,i),fn
        end do
      end do
    end do
  else
    write(50,'(4I6," : grid size, number of states")') np3d(:),ist1-ist0+1
    i=0
    do i3=0,ngridk(3)-1
      do i2=0,ngridk(2)-1
        do i1=0,ngridk(1)-1
          i=i+1
          ik=ivkik(i1,i2,i3)
          write(50,'(3G18.10)',advance='NO') vpc(:,i)
          do ist=ist0,ist1
            if (task == 101) then
! write the eigenvalues minus the Fermi energy separately
              write(50,'(F14.8)',advance='NO') evalsv(ist,ik)-efermi
            else
! write separate smooth delta functions at the Fermi energy
              x=(evalsv(ist,ik)-efermi)/swidth
              write(50,'(F14.8)',advance='NO') sdelta(stype,x)/swidth
            end if
          end do
          write(50,*)
        end do
      end do
    end do
  end if
  close(50)
end do
write(*,*)
write(*,'("Info(fermisurf):")')
if (ndmag == 1) then
  write(*,'(" 3D Fermi surface data written to FERMISURF_UP.OUT and &
   &FERMISURF_DN.OUT")')
else
  write(*,'(" 3D Fermi surface data written to FERMISURF.OUT")')
end if
if (task == 100) then
  write(*,'(" in terms of the product of eigenvalues minus the Fermi energy")')
else if (task == 101) then
  write(*,'(" in terms of separate eigenvalues minus the Fermi energy")')
else if (task == 103) then
  write(*,'(" in terms of a smooth delta function at the Fermi energy")')
else
  write(*,'(" in terms of separate delta functions at the Fermi energy")')
end if
deallocate(vpc)
! restore original parameters
reducek=reducek0
end subroutine

