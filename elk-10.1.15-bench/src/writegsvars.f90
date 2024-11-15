
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writegsvars
use modmain
use modvars
implicit none
call writevars('engytot',rv=engytot)
call writevars('engyx',rv=engyx)
call writevars('engyc',rv=engyc)
call writevars('efermi',rv=efermi)
call writevars('evalsum',rv=evalsum)
call writevars('engykn',rv=engykn)
call writevars('fermidos',rv=fermidos)
call writevars('bandgap',nv=2,rva=bandgap)
if (spinpol) then
  call writevars('momtot',nv=ndmag,rva=momtot)
  call writevars('momtotm',rv=momtotm)
  call writevars('mommt',nv=3*natmtot,rva=mommt)
end if
if (tforce) then
  call writevars('forcetot',nv=3*natmtot,rva=forcetot)
end if
end subroutine

