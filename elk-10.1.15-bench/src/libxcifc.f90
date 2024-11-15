
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module libxcifc

use xc_f03_lib_m

! libxc version number
integer libxcv(3)

private grad

contains

!BOP
! !ROUTINE: xcifc_libxc
! !INTERFACE:
subroutine xcifc_libxc(xctype,n,tempa,rho,rhoup,rhodn,g2rho,g2up,g2dn,grho2, &
 gup2,gdn2,gupdn,tau,tauup,taudn,ex,ec,vx,vc,vxup,vxdn,vcup,vcdn,dxdgr2,dxdgu2,&
 dxdgd2,dxdgud,dcdgr2,dcdgu2,dcdgd2,dcdgud,dxdg2r,dxdg2u,dxdg2d,dcdg2r,dcdg2u, &
 dcdg2d,wx,wxup,wxdn,wc,wcup,wcdn)
! !USES:
use, intrinsic :: iso_c_binding
! !INPUT/OUTPUT PARAMETERS:
!   xctype : type of exchange-correlation functional (in,integer(3))
!   n      : number of density points (in,integer)
!   tempa  : temperature in atomic units (in,real,optional)
!   rho    : spin-unpolarised charge density (in,real(n),optional)
!   rhoup  : spin-up charge density (in,real(n),optional)
!   rhodn  : spin-down charge density (in,real(n),optional)
!   g2rho  : grad^2 rho (in,real(n),optional)
!   g2up   : grad^2 rhoup (in,real(n),optional)
!   g2dn   : grad^2 rhodn (in,real(n),optional)
!   grho2  : |grad rho|^2 (in,real(n),optional)
!   gup2   : |grad rhoup|^2 (in,real(n),optional)
!   gdn2   : |grad rhodn|^2 (in,real(n),optional)
!   gupdn  : (grad rhoup).(grad rhodn) (in,real(n),optional)
!   tau    : kinetic energy density (in,real(n),optional)
!   tauup  : spin-up kinetic energy density (in,real(n),optional)
!   taudn  : spin-down kinetic energy density (in,real(n),optional)
!   ex     : exchange energy density (out,real(n),optional)
!   ec     : correlation energy density (out,real(n),optional)
!   vx     : spin-unpolarised exchange potential (out,real(n),optional)
!   vc     : spin-unpolarised correlation potential (out,real(n),optional)
!   vxup   : spin-up exchange potential (out,real(n),optional)
!   vxdn   : spin-down exchange potential (out,real(n),optional)
!   vcup   : spin-up correlation potential (out,real(n),optional)
!   vcdn   : spin-down correlation potential (out,real(n),optional)
!   dxdgr2 : de_x/d(|grad rho|^2) (out,real(n),optional)
!   dxdgu2 : de_x/d(|grad rhoup|^2) (out,real(n),optional)
!   dxdgd2 : de_x/d(|grad rhodn|^2) (out,real(n),optional)
!   dxdgud : de_x/d((grad rhoup).(grad rhodn)) (out,real(n),optional)
!   dcdgr2 : de_c/d(|grad rho|^2) (out,real(n),optional)
!   dcdgu2 : de_c/d(|grad rhoup|^2) (out,real(n),optional)
!   dcdgd2 : de_c/d(|grad rhodn|^2) (out,real(n),optional)
!   dcdgud : de_c/d((grad rhoup).(grad rhodn)) (out,real(n),optional)
!   dxdg2r : de_x/d(grad^2 rho) (out,real(n),optional)
!   dxdg2u : de_x/d(grad^2 rhoup) (out,real(n),optional)
!   dxdg2d : de_x/d(grad^2 rhodn) (out,real(n),optional)
!   dcdg2r : de_c/d(grad^2 rho) (out,real(n),optional)
!   dcdg2u : de_c/d(grad^2 rhoup) (out,real(n),optional)
!   dcdg2d : de_c/d(grad^2 rhodn) (out,real(n),optional)
!   wx     : de_x/dtau (out,real(n),optional)
!   wxup   : de_x/dtauup (out,real(n),optional)
!   wxdn   : de_x/dtaudn (out,real(n),optional)
!   wc     : de_c/dtau (out,real(n),optional)
!   wcup   : de_c/dtauup (out,real(n),optional)
!   wcdn   : de_c/dtaudn (out,real(n),optional)
! !DESCRIPTION:
!   Interface to the Libxc exchange-correlation functional library. The second
!   and third integers in {\tt xctype} define the exchange and correlation
!   functionals, respectively.
!
! !REVISION HISTORY:
!   Created April 2009 (Tyrel McQueen)
!   Modified September 2009 (JKD and TMQ)
!   Updated for Libxc 1, July 2010 (JKD)
!   Updated for Libxc 4, March 2018 (JKD)
!   Updated for Libxc 5, May 2020 (JKD)
!   Updated for Libxc 6, December 2022 (JKD)
!   Updated for Fortran 2003 interface, October 2024 (JKD)
!EOP
!BOC
implicit none
! mandatory arguments
integer, intent(in) :: xctype(3),n
! optional arguments
real(8), optional, intent(in) :: tempa
real(8), optional, intent(in) :: rho(n),rhoup(n),rhodn(n)
real(8), optional, intent(in) :: g2rho(n),g2up(n),g2dn(n)
real(8), optional, intent(in) :: grho2(n),gup2(n),gdn2(n),gupdn(n)
real(8), optional, intent(in) :: tau(n),tauup(n),taudn(n)
real(8), optional, intent(out) :: ex(n),ec(n),vx(n),vc(n)
real(8), optional, intent(out) :: vxup(n),vxdn(n),vcup(n),vcdn(n)
real(8), optional, intent(out) :: dxdgr2(n),dxdgu2(n),dxdgd2(n),dxdgud(n)
real(8), optional, intent(out) :: dxdg2r(n),dxdg2u(n),dxdg2d(n)
real(8), optional, intent(out) :: wx(n),wxup(n),wxdn(n)
real(8), optional, intent(out) :: dcdgr2(n),dcdgu2(n),dcdgd2(n),dcdgud(n)
real(8), optional, intent(out) :: dcdg2r(n),dcdg2u(n),dcdg2d(n)
real(8), optional, intent(out) :: wc(n),wcup(n),wcdn(n)
! local variables
integer nspin,fmly,id,k
integer(c_size_t) np
type(xc_f03_func_t) p
! allocatable arrays
real(8), allocatable :: r(:,:),sigma(:,:),vrho(:,:),vsigma(:,:)
real(8), allocatable :: lapl(:,:),t(:,:),vlapl(:,:),vtau(:,:)
if (present(rho)) then
  nspin=XC_UNPOLARIZED
else if (present(rhoup).and.present(rhodn)) then
  nspin=XC_POLARIZED
else
  write(*,*)
  write(*,'("Error(libxcifc): missing arguments")')
  write(*,*)
  stop
end if
if (xctype(2) /= 0) then
  if (xctype(2) == xctype(3)) then
    write(*,*)
    write(*,'("Error(libxcifc): Libxc exchange and correlation functionals")')
    write(*,'(" are the same : ",2I8)') xctype(2:3)
    write(*,*)
    stop
  end if
end if
! convert number of points to long integer
np=n
! loop over functional kinds (exchange or correlation)
do k=2,3
  id=xctype(k)
  if (id > 0) then
    fmly=xc_f03_family_from_id(id)
! initialise functional
    call xc_f03_func_init(p,id,nspin)
    select case(fmly)
    case(XC_FAMILY_LDA)
!-------------------------!
!     LDA functionals     !
!-------------------------!
! set temperature for free energy functional
      if ((id == XC_LDA_XC_KSDT).or.(id == XC_LDA_XC_GDSMFB)) then
        call xc_f03_func_set_ext_params(p,[tempa])
      end if
      if (k == 2) then
! exchange or a kinetic energy functional
        if (present(rho)) then
          if (present(ex)) then
            call xc_f03_lda_exc_vxc(p,np,rho,ex,vx)
          else
            call xc_f03_lda_vxc(p,np,rho,vx)
          end if
        else
          allocate(r(2,n),vrho(2,n))
          r(1,1:n)=rhoup(1:n); r(2,1:n)=rhodn(1:n)
          if (present(ex)) then
            call xc_f03_lda_exc_vxc(p,np,r,ex,vrho)
          else
            call xc_f03_lda_vxc(p,np,r,vrho)
          end if
          vxup(1:n)=vrho(1,1:n); vxdn(1:n)=vrho(2,1:n)
          deallocate(r,vrho)
        end if
      else
! correlation
        if (present(rho)) then
          if (present(ec)) then
            call xc_f03_lda_exc_vxc(p,np,rho,ec,vc)
          else
            call xc_f03_lda_vxc(p,np,rho,vc)
          end if
        else
          allocate(r(2,n),vrho(2,n))
          r(1,1:n)=rhoup(1:n); r(2,1:n)=rhodn(1:n)
          if (present(ec)) then
            call xc_f03_lda_exc_vxc(p,np,r,ec,vrho)
          else
            call xc_f03_lda_vxc(p,np,r,vrho)
          end if
          vcup(1:n)=vrho(1,1:n); vcdn=vrho(2,1:n)
          deallocate(r,vrho)
        end if
      end if
    case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
!-------------------------!
!     GGA functionals     !
!-------------------------!
      if (k == 2) then
! exchange or a kinetic energy functional
        if (present(rho)) then
          if (present(ex)) then
            call xc_f03_gga_exc_vxc(p,np,rho,grho2,ex,vx,dxdgr2)
          else
            call xc_f03_gga_vxc(p,np,rho,grho2,vx,dxdgr2)
          end if
        else
          allocate(r(2,n),sigma(3,n),vrho(2,n),vsigma(3,n))
          r(1,1:n)=rhoup(1:n); r(2,1:n)=rhodn(1:n)
          sigma(1,1:n)=gup2(1:n)
          if (present(gupdn)) then
            sigma(2,1:n)=gupdn(1:n)
          else
            sigma(2,1:n)=0.d0
          end if
          sigma(3,1:n)=gdn2(1:n)
          if (present(ex)) then
            call xc_f03_gga_exc_vxc(p,np,r,sigma,ex,vrho,vsigma)
          else
            call xc_f03_gga_vxc(p,np,r,sigma,vrho,vsigma)
          end if
          vxup(1:n)=vrho(1,1:n); vxdn(1:n)=vrho(2,1:n)
          dxdgu2(1:n)=vsigma(1,1:n)
          if (present(dxdgud)) dxdgud(1:n)=vsigma(2,1:n)
          dxdgd2(1:n)=vsigma(3,1:n)
          deallocate(r,sigma,vrho,vsigma)
        end if
      else
! correlation
        if (present(rho)) then
          if (present(ec)) then
            call xc_f03_gga_exc_vxc(p,np,rho,grho2,ec,vc,dcdgr2)
          else
            call xc_f03_gga_vxc(p,np,rho,grho2,vc,dcdgr2)
          end if
        else
          allocate(r(2,n),sigma(3,n),vrho(2,n),vsigma(3,n))
          r(1,1:n)=rhoup(1:n); r(2,1:n)=rhodn(1:n)
          sigma(1,1:n)=gup2(1:n)
          sigma(2,1:n)=gupdn(1:n)
          sigma(3,1:n)=gdn2(1:n)
          if (present(ec)) then
            call xc_f03_gga_exc_vxc(p,np,r,sigma,ec,vrho,vsigma)
          else
            call xc_f03_gga_vxc(p,np,r,sigma,vrho,vsigma)
          end if
          vcup(1:n)=vrho(1,1:n); vcdn(1:n)=vrho(2,1:n)
          dcdgu2(1:n)=vsigma(1,1:n)
          dcdgud(1:n)=vsigma(2,1:n)
          dcdgd2(1:n)=vsigma(3,1:n)
          deallocate(r,sigma,vrho,vsigma)
        end if
      end if
    case(XC_FAMILY_MGGA)
!------------------------------!
!     meta-GGA functionals     !
!------------------------------!
      if (k == 2) then
! exchange or a kinetic energy functional
        if (present(rho)) then
          if (present(ex)) then
            call xc_f03_mgga_exc_vxc(p,np,rho,grho2,g2rho,tau,ex,vx,dxdgr2, &
             dxdg2r,wx)
          else
            call xc_f03_mgga_vxc(p,np,rho,grho2,g2rho,tau,vx,dxdgr2,dxdg2r,wx)
          end if
        else
          allocate(r(2,n),sigma(3,n),lapl(2,n),t(2,n))
          allocate(vrho(2,n),vsigma(3,n),vlapl(2,n),vtau(2,n))
          r(1,1:n)=rhoup(1:n); r(2,1:n)=rhodn(1:n)
          sigma(1,1:n)=gup2(1:n); sigma(3,1:n)=gdn2(1:n)
          if (present(gupdn)) then
            sigma(2,1:n)=gupdn(1:n)
          else
            sigma(2,1:n)=0.d0
          end if
          lapl(1,1:n)=g2up(1:n); lapl(2,1:n)=g2dn(1:n)
          t(1,1:n)=tauup(1:n); t(2,1:n)=taudn(1:n)
          if (present(ex)) then
            call xc_f03_mgga_exc_vxc(p,np,r,sigma,lapl,t,ex,vrho,vsigma,vlapl, &
             vtau)
          else
            call xc_f03_mgga_vxc(p,np,r,sigma,lapl,t,vrho,vsigma,vlapl,vtau)
          end if
          vxup(1:n)=vrho(1,1:n); vxdn(1:n)=vrho(2,1:n)
          dxdgu2(1:n)=vsigma(1,1:n); dxdgd2(1:n)=vsigma(3,1:n)
          if (present(dxdgud)) dxdgud(1:n)=vsigma(2,1:n)
          dxdg2u(1:n)=vlapl(1,1:n); dxdg2d(1:n)=vlapl(2,1:n)
          wxup(1:n)=vtau(1,1:n); wxdn(1:n)=vtau(2,1:n)
          deallocate(r,sigma,lapl,t)
          deallocate(vrho,vsigma,vlapl,vtau)
        end if
      else
! correlation
        if (present(rho)) then
          if (present(ec)) then
            call xc_f03_mgga_exc_vxc(p,np,rho,grho2,g2rho,tau,ec,vc,dcdgr2, &
             dcdg2r,wc)
          else
            call xc_f03_mgga_vxc(p,np,rho,grho2,g2rho,tau,vc,dcdgr2,dcdg2r,wc)
          end if
        else
          allocate(r(2,n),sigma(3,n),lapl(2,n),t(2,n))
          allocate(vrho(2,n),vsigma(3,n),vlapl(2,n),vtau(2,n))
          r(1,1:n)=rhoup(1:n); r(2,1:n)=rhodn(1:n)
          sigma(1,1:n)=gup2(1:n)
          sigma(2,1:n)=gupdn(1:n)
          sigma(3,1:n)=gdn2(1:n)
          lapl(1,1:n)=g2up(1:n); lapl(2,1:n)=g2dn(1:n)
          t(1,1:n)=tauup(1:n); t(2,1:n)=taudn(1:n)
          if (present(ec)) then
            call xc_f03_mgga_exc_vxc(p,np,r,sigma,lapl,t,ec,vrho,vsigma,vlapl, &
             vtau)
          else
            call xc_f03_mgga_vxc(p,np,r,sigma,lapl,t,vrho,vsigma,vlapl,vtau)
          end if
          vcup(1:n)=vrho(1,1:n); vcdn(1:n)=vrho(2,1:n)
          dcdgu2(1:n)=vsigma(1,1:n)
          dcdgud(1:n)=vsigma(2,1:n)
          dcdgd2(1:n)=vsigma(3,1:n)
          dcdg2u(1:n)=vlapl(1,1:n); dcdg2d(1:n)=vlapl(2,1:n)
          wcup(1:n)=vtau(1,1:n); wcdn(1:n)=vtau(2,1:n)
          deallocate(r,sigma,lapl,t)
          deallocate(vrho,vsigma,vlapl,vtau)
        end if
      end if
    case default
      write(*,*)
      write(*,'("Error(libxcifc): unsupported Libxc functional family : ",I8)')&
       fmly
      write(*,*)
      stop
    end select
! destroy functional
    call xc_f03_func_end(p)
  else
! case when id=0
    if (k == 2) then
      if (present(ex)) ex(1:n)=0.d0
      if (present(vx)) vx(1:n)=0.d0
      if (present(vxup)) vxup(1:n)=0.d0
      if (present(vxdn)) vxdn(1:n)=0.d0
      if (present(dxdgr2)) dxdgr2(1:n)=0.d0
      if (present(dxdgu2)) dxdgu2(1:n)=0.d0
      if (present(dxdgd2)) dxdgd2(1:n)=0.d0
      if (present(dxdgud)) dxdgud(1:n)=0.d0
      if (present(dxdg2r)) dxdg2r(1:n)=0.d0
      if (present(dxdg2u)) dxdg2u(1:n)=0.d0
      if (present(dxdg2d)) dxdg2d(1:n)=0.d0
      if (present(wx)) wx(1:n)=0.d0
      if (present(wxup)) wxup(1:n)=0.d0
      if (present(wxdn)) wxdn(1:n)=0.d0
    else
      if (present(ec)) ec(1:n)=0.d0
      if (present(vc)) vc(1:n)=0.d0
      if (present(vcup)) vcup(1:n)=0.d0
      if (present(vcdn)) vcdn(1:n)=0.d0
      if (present(dcdgr2)) dcdgr2(1:n)=0.d0
      if (present(dcdgu2)) dcdgu2(1:n)=0.d0
      if (present(dcdgd2)) dcdgd2(1:n)=0.d0
      if (present(dcdgud)) dcdgud(1:n)=0.d0
      if (present(dcdg2r)) dcdg2r(1:n)=0.d0
      if (present(dcdg2u)) dcdg2u(1:n)=0.d0
      if (present(dcdg2d)) dcdg2d(1:n)=0.d0
      if (present(wc)) wc(1:n)=0.d0
      if (present(wcup)) wcup(1:n)=0.d0
      if (present(wcdn)) wcdn(1:n)=0.d0
    end if
  end if
end do
end subroutine

subroutine fxcifc_libxc(fxctype,n,rho,rhoup,rhodn,fxc,fxcuu,fxcud,fxcdd)
use, intrinsic :: iso_c_binding
implicit none
! mandatory arguments
integer, intent(in) :: fxctype(3),n
! optional arguments
real(8), optional, intent(in) :: rho(n),rhoup(n),rhodn(n)
real(8), optional, intent(out) :: fxc(n),fxcuu(n),fxcud(n),fxcdd(n)
! local variables
integer nspin,fmly,id,k
integer(c_size_t) np
type(xc_f03_func_t) p
! allocatable arrays
real(8), allocatable :: r(:,:),f(:,:)
np=n
if (present(rho)) then
  nspin=XC_UNPOLARIZED
else if (present(rhoup).and.present(rhodn)) then
  nspin=XC_POLARIZED
else
  write(*,*)
  write(*,'("Error(libxcifc): missing arguments")')
  write(*,*)
  stop
end if
! zero the kernel
if (present(fxc)) fxc(1:n)=0.d0
if (present(fxcuu)) fxcuu(1:n)=0.d0
if (present(fxcud)) fxcud(1:n)=0.d0
if (present(fxcdd)) fxcdd(1:n)=0.d0
! loop over functional kinds (exchange or correlation)
do k=2,3
  id=fxctype(k)
  if (id <= 0) cycle
  fmly=xc_f03_family_from_id(id)
! initialise functional
  call xc_f03_func_init(p,id,nspin)
  select case(fmly)
  case(XC_FAMILY_LDA)
!-------------------------!
!     LDA functionals     !
!-------------------------!
    if (present(rho)) then
      allocate(f(1,n))
      call xc_f03_lda_fxc(p,np,rho,f)
      fxc(1:n)=fxc(1:n)+f(1,1:n)
      deallocate(f)
    else
      allocate(r(2,n),f(3,n))
      r(1,1:n)=rhoup(1:n); r(2,1:n)=rhodn(1:n)
      call xc_f03_lda_fxc(p,np,r,f)
      fxcuu(1:n)=fxcuu(1:n)+f(1,1:n)
      fxcud(1:n)=fxcud(1:n)+f(2,1:n)
      fxcdd(1:n)=fxcdd(1:n)+f(3,1:n)
      deallocate(r,f)
    end if
  case default
    write(*,*)
    write(*,'("Error(libxcifc): unsupported Libxc functional family : ",I8)') &
     fmly
    write(*,'(" for calculating f_xc")')
    write(*,*)
    stop
  end select
! destroy functional
  call xc_f03_func_end(p)
end do
end subroutine

subroutine xcdata_libxc(xctype,xcdescr,xcspin,xcgrad,hybrid,hybridc)
implicit none
! arguments
integer, intent(in) :: xctype(3)
character(264), intent(out) :: xcdescr
integer, intent(out) :: xcspin,xcgrad
logical, intent(out) :: hybrid
real(8), intent(out) :: hybridc
! local variables
integer k,id,fmly,g
character(128) name
type(xc_f03_func_t) p
! check version is compatible
call xc_f03_version(libxcv(1),libxcv(2),libxcv(3))
if (all(libxcv(1) /= [6,7])) then
  write(*,*)
  write(*,'("Error(libxcifc): incompatible Libxc version : ",I0,".",I0,".",&
   &I0)') libxcv(:)
  write(*,*)
  stop
end if
! undefined spin polarisation
xcspin=-1
! undefined gradient type
xcgrad=-1
! not a hybrid functional by default
hybrid=.false.
! default description
xcdescr='none'
do k=2,3
  id=xctype(k)
  if (id <= 0) cycle
! initialise functional
  call xc_f03_func_init(p,id,XC_UNPOLARIZED)
!-----------------------!
!     gradient type     !
!-----------------------!
  g=grad(p)
  if (g < 0) then
    write(*,*)
    write(*,'("Error(libxcifc): unsupported gradient type")')
    write(*,*)
    stop
  end if
  if (k == 3) then
    if ((xcgrad >= 0).and.(g /= xcgrad)) then
      write(*,*)
      write(*,'("Error(libxcifc): inconsistent exchange and correlation &
       &gradient types")')
      write(*,*)
      stop
    end if
  end if
  xcgrad=g
!----------------------------!
!     hybrid functionals     !
!----------------------------!
  fmly=xc_f03_family_from_id(id)
! hybrid GGA functionals
  if (fmly == XC_FAMILY_HYB_GGA) then
    if (k == 2) then
      write(*,*)
      write(*,'("Error(libxcifc): set only correlation part of xctype for &
       &Libxc hybrids")')
      write(*,*)
      stop
    end if
! set the hybrid mixing coefficient
    hybridc=xc_f03_hyb_exx_coef(p)
    hybrid=.true.
  end if
!--------------------------------!
!     functional description     !
!--------------------------------!
  name=xc_f03_functional_get_name(id)
  if (k == 2) then
    xcdescr=trim(name)
  else
    xcdescr=trim(xcdescr)//'; '//trim(name)
  end if
! destroy functional
  call xc_f03_func_end(p)
end do
end subroutine

integer function grad(p)
implicit none
! arguments
type(xc_f03_func_t), intent(in) :: p
! local variables
type(xc_f03_func_info_t) info
type(xc_f03_func_t) pa
integer fmly,flg,knd
integer naux,i
! allocatable arrays
integer, allocatable :: ida(:)
info=xc_f03_func_get_info(p)
fmly=xc_f03_func_info_get_family(info)
select case(fmly)
case(XC_FAMILY_LDA)
! no gradients required
  grad=0
case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
! GGA gradients required
  grad=2
case(XC_FAMILY_MGGA)
  flg=xc_f03_func_info_get_flags(info)
  if (iand(flg,XC_FLAGS_NEEDS_LAPLACIAN) == 0) then
! tau and no Laplacian
    grad=4
  else
! tau and Laplacian
    grad=5
  end if
! check if tau is not required
  knd=xc_f03_func_info_get_kind(info)
  if (knd == XC_KINETIC) then
! this is a kinetic energy functional so tau is not required
    grad=3
  else
! this is an exchange and/or correlation functional
    naux=xc_f03_num_aux_funcs(p)
    allocate(ida(naux))
    call xc_f03_aux_func_ids(p,ida)
    do i=1,naux
      call xc_f03_func_init(pa,ida(i),XC_UNPOLARIZED)
      info=xc_f03_func_get_info(pa)
      knd=xc_f03_func_info_get_kind(info)
      call xc_f03_func_end(pa)
      if (knd == XC_KINETIC) then
! this is a deorbitalised functional and does not require tau
        grad=3
        exit
      end if
    end do
    deallocate(ida)
  end if
case default
  grad=-1
end select
end function
!EOC

end module

