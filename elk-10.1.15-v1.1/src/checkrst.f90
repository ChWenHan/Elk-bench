
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine checkrst
use modmain
use modmpi
use moddelf
implicit none
! check for RESTART file (only MPI master process)
if (mp_mpi) then
  inquire(file='RESTART',exist=trestart)
  if (trestart) then
    write(*,'("Info(checkrst): RESTART file exists")')
! delete the RESTART file
    call delfile('RESTART')
  end if
end if
! broadcast trestart from master process to all other processes
call mpi_bcast(trestart,1,mpi_logical,0,mpicom,ierror)
end subroutine

