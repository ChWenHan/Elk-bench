
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writephn
use modmain
use modphonon
implicit none
! local variables
integer iq,i,j,is,ia,ip
! allocatable arrays
real(8), allocatable :: wq(:)
complex(8), allocatable :: dq(:,:),ev(:,:)
! initialise universal variables
call init0
call init2
call initph
allocate(wq(nbph),dq(nbph,nbph),ev(nbph,nbph))
open(50,file='PHONON.OUT',form='FORMATTED')
do iq=1,nphwrt
  call dynrtoq(vqlwrt(:,iq),dynr,dq)
  call dynev(dq,wq,ev)
  write(50,*)
  write(50,'(I6,3G18.10," : q-point, vqlwrt")') iq,vqlwrt(:,iq)
  do j=1,nbph
    write(50,*)
    write(50,'(I6,G18.10," : mode, frequency")') j,wq(j)
    i=0
    do is=1,nspecies
      do ia=1,natoms(is)
        do ip=1,3
          i=i+1
          if (i == 1) then
            write(50,'(3I4,2G18.10," : species, atom, polarisation, &
             &eigenvector")') is,ia,ip,ev(i,j)
          else
            write(50,'(3I4,2G18.10)') is,ia,ip,ev(i,j)
          end if
        end do
      end do
    end do
  end do
  write(50,*)
end do
close(50)
write(*,*)
write(*,'("Info(writephn): phonon frequencies and eigenvectors written to &
 &PHONON.OUT")')
write(*,'(" for all q-vectors in the phwrite list")')
deallocate(wq,dq,ev)
end subroutine

