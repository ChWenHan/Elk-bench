
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeepsinv
use modmain
use modmpi
implicit none
! initialise global variables
call init0
call init1
call init2
call init3
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! read the eigenvalues and occupation numbers from file
call readevalsv
call readoccsv
! generate the inverse dielectric function and write to file
call epsinv
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(writeepsinv):")')
  write(*,'(" inverse RPA dielectric function, ϵ⁻¹(G,G'',q,w), written to &
   &EPSINV.OUT")')
end if
end subroutine

