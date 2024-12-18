
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rdmengyxc
! !INTERFACE:
subroutine rdmengyxc
! !USES:
use modmain
use modrdm
! !DESCRIPTION:
!   Calculates RDMFT exchange-correlation energy.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik1,ik2,jk,iv(3)
integer ist1,ist2
real(8) t1,t2,t3,t4
! allocatable arays
real(8), allocatable :: vcl1221(:,:,:)
! calculate the prefactor
if (rdmxctype == 0) then
  engyx=0.d0
  return
else if (rdmxctype == 1) then
! Hartree-Fock functional
  t1=0.5d0/occmax
else if (rdmxctype == 2) then
! Power functional
  if (spinpol) then
    t1=0.5d0
  else
    t1=(0.25d0)**rdmalpha
  end if
else
  write(*,*)
  write(*,'("Error(rdmengyxc): rdmxctype not defined : ",I8)') rdmxctype
  write(*,*)
  stop
end if
! exchange-correlation energy
engyx=0.d0
allocate(vcl1221(nstsv,nstsv,nkpt))
! start loop over non-reduced k-points
do ik1=1,nkptnr
  call getvcl1221(ik1,vcl1221)
! find the equivalent reduced k-point
  iv(:)=ivk(:,ik1)
  jk=ivkik(iv(1),iv(2),iv(3))
  do ist1=1,nstsv
! start loop over reduced k-points
   do ik2=1,nkpt
     do ist2=1,nstsv
! Hartree-Fock functional
        if (rdmxctype == 1) then
          t2=t1*wkpt(ik2)*occsv(ist2,ik2)*occsv(ist1,jk)
! Power functional
        else if (rdmxctype == 2) then
          t3=occsv(ist2,ik2)*occsv(ist1,jk)
          t4=sum(abs(vkl(:,ik2)-vkl(:,ik1)))
          if ((ist2 == ist1).and.(t4 < epslat)) then
            t2=(0.5d0/occmax)*wkpt(ik2)*t3
          else
            t2=t1*wkpt(ik2)*(t3**rdmalpha)
          end if
        end if
        engyx=engyx-t2*vcl1221(ist2,ist1,ik2)
      end do
    end do
  end do
end do
deallocate(vcl1221)
end subroutine
!EOC

