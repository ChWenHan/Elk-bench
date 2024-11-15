
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tddftlr
use modmain
use modtddft
use modtest
use modmpi
use modomp
implicit none
! local variables
logical tq0
integer, parameter :: maxit=500
integer iq,ik,isym
integer nm,it,i,j,n
integer iw,ioc,nthd
real(8) v(3),t1,t2
complex(8) vfxcp,z1
character(256) fname
! allocatable arrays
integer(omp_lock_kind), allocatable :: lock(:)
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: vchi0(:,:,:),vfxc(:,:,:)
complex(8), allocatable :: eps0(:,:,:),epsi(:,:,:),epsm(:,:,:)
complex(8), allocatable :: zw(:),a(:,:)
! initialise global variables
call init0
call init1
call init2
call init3
! check q-vector is commensurate with k-point grid
v(1:3)=dble(ngridk(1:3))*vecql(1:3)
v(1:3)=abs(v(1:3)-nint(v(1:3)))
if ((v(1) > epslat).or.(v(2) > epslat).or.(v(3) > epslat)) then
  write(*,*)
  write(*,'("Error(tddftlr): q-vector incommensurate with k-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" vecql : ",3G18.10)') vecql
  write(*,*)
  stop
end if
! find the equivalent reduced q-point
call findqpt(vecql,isym,iq)
! check if q = 0
tq0=.false.
if (sum(abs(vecql(:))) < epslat) tq0=.true.
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the eigenvalues and occupation numbers from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
! generate the G+q-vectors and related functions
allocate(vgqc(3,ngrf),gqc(ngrf),jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
call gengqf(ngrf,vecqc,vgqc,gqc,jlgqr,ylmgq,sfacgq)
deallocate(vgqc)
! generate the regularised Coulomb Green's function in G+q-space
allocate(gclgq(ngrf))
call gengclgq(.true.,iq,ngrf,gqc,gclgq)
gclgq(1:ngrf)=sqrt(gclgq(1:ngrf))
! matrix sizes
if (tq0) then
! for q = 0 the head is a 3 x 3 matrix and the wings are 3 x ngrf
  nm=ngrf+2
else
! otherwise the head is just G = G' = 0 and finite q
  nm=ngrf
end if
! allocate local arrays
allocate(vchi0(nm,nm,nwrf),vfxc(nm,nm,nwrf))
allocate(eps0(nm,nm,nwrf),epsi(nm,nm,nwrf))
! initialise the OpenMP locks
allocate(lock(nwrf))
do iw=1,nwrf
  call omp_init_lock(lock(iw))
end do
! compute v¹⸍² χ₀ v¹⸍² (the symmetric version of v χ₀)
vchi0(1:nm,1:nm,1:nwrf)=0.d0
call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
!$OMP CRITICAL(tddftlr_)
  write(*,'("Info(tddftlr): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL(tddftlr_)
! add to v¹⸍² χ₀ v¹⸍²
  call genvchi0(.true.,ik,lock,vecql,gclgq,jlgqr,ylmgq,sfacgq,nm,vchi0)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! destroy the OpenMP locks
do iw=1,nwrf
  call omp_destroy_lock(lock(iw))
end do
deallocate(lock)
! add vchi0 from each process and redistribute
if (np_mpi > 1) then
  n=nm*nm*nwrf
  call mpi_allreduce(mpi_in_place,vchi0,n,mpi_double_complex,mpi_sum,mpicom, &
   ierror)
end if
! calculate symmetric ϵ_0 = 1 - v¹⸍² χ₀ v¹⸍²
eps0(1:nm,1:nm,1:nwrf)=-vchi0(1:nm,1:nm,1:nwrf)
do i=1,nm
  eps0(i,i,1:nwrf)=eps0(i,i,1:nwrf)+1.d0
end do
! initialise ϵ for use with the bootstrap functional
if (any(fxctype(1) == [210,211])) then
  epsi(1:nm,1:nm,1:nwrf)=vchi0(1:nm,1:nm,1:nwrf)
  do i=1,nm
    epsi(i,i,1:nwrf)=epsi(i,i,1:nwrf)+1.d0
  end do
end if
allocate(a(nm,nm))
vfxcp=0.d0
it=0
10 continue
! compute v χ₀ v⁻¹⸍² f_xc v⁻¹⸍² v χ₀
call genvfxc(tq0,.true.,gclgq,nm,vchi0,eps0,epsi,vfxc)
! begin loop over frequencies
do iw=1,nwrf
! compute 1 - v¹⸍² χ₀ v¹⸍² - v⁻¹⸍² f_xc v⁻¹⸍² v χ₀
  a(1:nm,1:nm)=eps0(1:nm,1:nm,iw)-vfxc(1:nm,1:nm,iw)
! invert this matrix
  call zminv(nm,a)
! left multiply by v¹⸍² χ₀ v¹⸍²
  call zgemm('N','N',nm,nm,nm,zone,vchi0(:,:,iw),nm,a,nm,zzero,epsi(:,:,iw),nm)
! compute ϵ⁻¹ = 1 + v¹⸍² χ v¹⸍²
  do i=1,nm
    epsi(i,i,iw)=1.d0+epsi(i,i,iw)
  end do
end do
if (fxctype(1) == 210) then
! self-consistent bootstrap f_xc
  it=it+1
  if (it > maxit) then
    write(*,*)
    write(*,'("Error(tddftlr): bootstrap kernel failed to converge")')
    write(*,*)
    stop
  end if
  if (mp_mpi) then
    if (mod(it,10) == 0) then
      write(*,'("Info(tddftlr): done ",I4," bootstrap iterations")') it
      t1=dble(vfxc(1,1,1))
      write(*,'(" head of matrix v⁻¹⸍² f_xc v⁻¹⸍² : ",G18.10)') t1
      write(*,'("  multiplied by -4π gives α : ",G18.10)') -fourpi*t1
    end if
  end if
! check for convergence
  t1=abs(vfxcp)-abs(vfxc(1,1,1))
  vfxcp=vfxc(1,1,1)
  if (abs(t1) > 1.d-8) goto 10
else if (fxctype(1) == 211) then
! single iteration bootstrap
  it=it+1
  if (it <= 1) goto 10
end if
deallocate(gclgq,jlgqr)
deallocate(ylmgq,sfacgq)
deallocate(vchi0,vfxc)
! invert ϵ⁻¹ to find ϵ and store in array eps0
do iw=1,nwrf
  eps0(1:nm,1:nm,iw)=epsi(1:nm,1:nm,iw)
  call zminv(nm,eps0(:,:,iw))
end do
if (mp_mpi) then
! write G = G' = 0 components to file
  if (tq0) then
    do ioc=1,noptcomp
      i=optcomp(1,ioc)
      j=optcomp(2,ioc)
      write(fname,'("EPSILON_TDDFT_",2I1,".OUT")') i,j
      open(50,file=trim(fname),form='FORMATTED')
      write(fname,'("EPSINV_TDDFT_",2I1,".OUT")') i,j
      open(51,file=trim(fname),form='FORMATTED')
      do iw=2,nwrf
        write(50,'(2G18.10)') dble(wrf(iw)),dble(eps0(i,j,iw))
        write(51,'(2G18.10)') dble(wrf(iw)),dble(epsi(i,j,iw))
      end do
      write(50,*)
      write(51,*)
      do iw=2,nwrf
        write(50,'(2G18.10)') dble(wrf(iw)),aimag(eps0(i,j,iw))
        write(51,'(2G18.10)') dble(wrf(iw)),aimag(epsi(i,j,iw))
      end do
      close(50)
      close(51)
    end do
  else
    open(50,file='EPSILON_TDDFT.OUT',form='FORMATTED')
    open(51,file='EPSINV_TDDFT.OUT',form='FORMATTED')
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),dble(eps0(1,1,iw))
      write(51,'(2G18.10)') dble(wrf(iw)),dble(epsi(1,1,iw))
    end do
    write(50,*)
    write(51,*)
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),aimag(eps0(1,1,iw))
      write(51,'(2G18.10)') dble(wrf(iw)),aimag(epsi(1,1,iw))
    end do
    close(50)
    close(51)
  end if
! find the macroscopic part of ϵ by inverting the 3x3 head only
  if (tq0) then
    allocate(epsm(3,3,nwrf))
    do iw=1,nwrf
      epsm(1:3,1:3,iw)=epsi(1:3,1:3,iw)
      call zminv(3,epsm(:,:,iw))
    end do
! write out the macroscopic components
    do ioc=1,noptcomp
      i=optcomp(1,ioc)
      j=optcomp(2,ioc)
      write(fname,'("EPSM_TDDFT_",2I1,".OUT")') i,j
      open(50,file=trim(fname),form='FORMATTED')
      do iw=2,nwrf
        write(50,'(2G18.10)') dble(wrf(iw)),dble(epsm(i,j,iw))
      end do
      write(50,*)
      do iw=2,nwrf
        write(50,'(2G18.10)') dble(wrf(iw)),aimag(epsm(i,j,iw))
      end do
      close(50)
    end do
    allocate(zw(nwrf))
! output the Faraday angle parameters Δδ and Δβ
    do iw=2,nwrf
      zw(iw)=0.5d0*zi*epsm(1,2,iw)/sqrt(epsm(1,1,iw))
    end do
    open(50,file='FARADAY.OUT',form='FORMATTED')
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),dble(zw(iw))
    end do
    write(50,*)
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),aimag(zw(iw))
    end do
    close(50)
! output the Kerr angle
    do iw=2,nwrf
      zw(iw)=-epsm(1,2,iw)/(sqrt(epsm(1,1,iw))*(epsm(1,1,iw)-1.d0))
    end do
    open(50,file='KERR_TDDFT.OUT',form='FORMATTED')
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),dble(zw(iw))*180.d0/pi
    end do
    write(50,*)
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),aimag(zw(iw))*180.d0/pi
    end do
    close(50)
! output magnetic linear dichroism (MLD) spectrum
    t1=sin(thetamld)**2
    t2=sin(2.d0*thetamld)
    do iw=2,nwrf
      z1=epsm(1,1,iw)
      zw(iw)=t2*epsm(1,2,iw)/((z1-1.d0)*(z1-(t1*(z1+1.d0))))
    end do
    open(50,file='MLD.OUT',form='FORMATTED')
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),dble(zw(iw))
    end do
    write(50,*)
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),aimag(zw(iw))
    end do
    close(50)
    deallocate(epsm,zw)
  end if
end if
! write inverse ϵ to test file
call writetest(320,'inverse ϵ',nv=nm*nm*nwrf,tol=1.d-2,zva=epsi)
deallocate(eps0,epsi,a)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(tddftlr):")')
  if (tq0) then
    write(*,'(" Dielectric tensor written to EPSILON_TDDFT_ij.OUT")')
    write(*,'(" Inverse written to EPSINV_TDDFT_ij.OUT")')
    write(*,'(" Macroscopic part written to EPSM_TDDFT_ij.OUT")')
    write(*,'(" for components")')
    do ioc=1,noptcomp
      write(*,'("  i = ",I1,", j = ",I1)') optcomp(1:2,ioc)
    end do
    write(*,*)
    write(*,'(" Faraday angle parameters Δδ and Δβ written to FARADAY.OUT")')
    write(*,'(" MOKE Kerr angle written to KERR_TDDFT.OUT")')
    write(*,'(" Magnetic linear dichroism (MLD) spectrum written to MLD.OUT")')
    write(*,*)
    write(*,'(" Note that the q-vector is zero and therefore the head of the")')
    write(*,'(" tensor is a 3 x 3 matrix and the wings are 3 x ngrf matrices")')
  else
    write(*,'(" Dielectric tensor written to EPSILON_TDDFT.OUT")')
    write(*,'(" Inverse written to EPSINV_TDDFT.OUT")')
    write(*,'(" q-vector (lattice coordinates) : ")')
    write(*,'(3G18.10)') vecql
    write(*,'(" q-vector length : ",G18.10)') gqc(1)
  end if
end if
deallocate(gqc)
end subroutine

