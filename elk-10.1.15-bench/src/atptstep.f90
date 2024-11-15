
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine atptstep(ft)
use modmain
use modtddft
use modmpi
implicit none
! arguments
real(8), intent(in) :: ft(3,natmtot)
! local variables
integer itn,is,ia,ias
real(8) dt,t1,t2
! next time step when the forces will be calculated
itn=itimes+ntsforce
if (itn > ntimes) return
! time step length
dt=times(itn)-times(itimes)
do is=1,nspecies
  t1=1.d0-atdfc*dt
  if (t1 < 0.d0) t1=0.d0
  t2=dt/spmass(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! add to the atomic velocities
    atdvc(1:3,1,ia,is)=t1*atdvc(1:3,1,ia,is)+t2*ft(1:3,ias)
! add to the atomic displacements
    atdvc(1:3,0,ia,is)=atdvc(1:3,0,ia,is)+atdvc(1:3,1,ia,is)*dt
  end do
end do
! write the atomic displacements and velocities to file
if (mp_mpi) call writeatdvc
end subroutine

