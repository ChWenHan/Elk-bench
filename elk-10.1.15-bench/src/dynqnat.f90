
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynqnat(sgn,vpl,dq)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: sgn
real(8), intent(in) :: vpl(3)
complex(8), intent(inout) :: dq(nbph,nbph)
! local variables
integer ias,jas,ip,jp,i,j
real(8) vpfl(3),vpc(3),v(3),t0,t1
! automatic arrays
real(8) vbpc(3,natmtot)
! map vpl to the first Brillouin zone
vpfl(:)=vpl(:)
call vecfbz(epslat,bvec,vpfl)
! q-vector in Cartesian coordinates
call r3mv(bvec,vpfl,vpc)
! multiply the q-vector with the Born effective charge tensor
do ias=1,natmtot
  call r3mv(bec(:,:,ias),vpc,vbpc(:,ias))
end do
! compute q⋅ϵ⋅q
call r3mv(epsw0,vpc,v)
t0=vpc(1)*v(1)+vpc(2)*v(2)+vpc(3)*v(3)
if (abs(t0) > 1.d-8) then
  t0=fourpi/(omega*t0)
  if (sgn < 0) t0=-t0
else
  t0=0.d0
end if
! ensure that the function tends to zero smoothly near the zone boundary
do i=1,3
  t1=abs(vpfl(i))
  if (t1 < 0.5d0) then
    t0=t0*(1.d0-2.d0*t1)
  else
    t0=0.d0
  end if
end do
j=0
do jas=1,natmtot
  do jp=1,3
    j=j+1
    i=0
    do ias=1,natmtot
      do ip=1,3
        i=i+1
        t1=vbpc(ip,ias)*vbpc(jp,jas)
        dq(i,j)=dq(i,j)+t0*t1
      end do
    end do
  end do
end do
end subroutine

