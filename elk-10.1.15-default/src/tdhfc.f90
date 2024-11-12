
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tdhfc
use modmain
use modtdhfc
implicit none
! initialise global variables
call init0
call init1
call readstate
call genvsig
call gencore
call linengy
call genapwlofr
call gensocfr
! read the eigenvalues and occupation numbers from file
call readevalsv
call readoccsv
! find the number of and index to the states within the energy window
call genidxthc


end subroutine

