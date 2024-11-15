
! Copyright (C) 2023 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine cpotcoul(nr,nri,np,ld1,rl,ngdg,igf,ngp,gpc,gclgp,ld2,jlgprmt,ylmgp, &
 sfacgp,crhoir,ld3,cvclmt,cvclir)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
integer, intent(in) :: ld1
real(8), intent(in) :: rl(ld1,-lmaxo-1:lmaxo+2,nspecies)
integer, intent(in) :: ngdg(3),igf(*),ngp
real(8), intent(in) :: gpc(ngp),gclgp(ngp)
integer, intent(in) :: ld2
real(8), intent(in) :: jlgprmt(0:lnpsd,ld2,nspecies)
complex(8), intent(in) :: ylmgp(lmmaxo,ngp),sfacgp(ld2,natmtot)
complex(4), intent(in) :: crhoir(*)
integer, intent(in) :: ld3
complex(4), intent(inout) :: cvclmt(ld3,natmtot)
complex(4), intent(out) :: cvclir(*)
! local variables
integer is,ia,ias,iro
integer l,lm,lma,lmb
integer ig,jg,i,i0,i1
real(8) t1,t2,t3
complex(8) z1,z2
! automatic arrays
real(8) rl2(0:lmaxo)
complex(8) qlm(lmmaxo,natmtot),zlm(lmmaxo)
! external functions
real(8), external :: factn2
! compute the multipole moments from the muffin-tin potentials
t1=1.d0/fourpi
do ias=1,natmtot
  is=idxis(ias)
  i=np(is)-lmmaxo
  do l=0,lmaxo
    t2=t1*dble(2*l+1)*rmtl(l+1,is)
    lma=l**2+1; lmb=lma+2*l
    qlm(lma:lmb,ias)=t2*cvclmt(i+lma:i+lmb,ias)
  end do
end do
! Fourier transform density to G-space and store in zvclir
call ccopy(ngdg(1)*ngdg(2)*ngdg(3),crhoir,1,cvclir,1)
call cfftifc(3,ngdg,-1,cvclir)
! subtract the multipole moments of the interstitial charge density
do is=1,nspecies
  rl2(0:lmaxo)=rmtl(2:lmaxo+2,is)
  t1=rl2(0)*fourpi*y00
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    zlm(1:lmmaxo)=0.d0
    do ig=1,ngp
      jg=igf(ig)
      if (gpc(ig) > epslat) then
        z1=cvclir(jg)*sfacgp(ig,ias)/gpc(ig)
        zlm(1)=zlm(1)+jlgprmt(1,ig,is)*t1*z1
        do l=1,lmaxo
          lma=l**2+1; lmb=lma+2*l
          z2=jlgprmt(l+1,ig,is)*rl2(l)*z1
          zlm(lma:lmb)=zlm(lma:lmb)+z2*conjg(ylmgp(lma:lmb,ig))
        end do
      else
        t2=(fourpi/3.d0)*rmtl(3,is)*y00
        zlm(1)=zlm(1)+t2*cvclir(jg)
      end if
    end do
    qlm(1:lmmaxo,ias)=qlm(1:lmmaxo,ias)-zlm(1:lmmaxo)
  end do
end do
! find the smooth pseudocharge within the muffin-tin whose multipoles are the
! difference between the real muffin-tin and interstitial multipoles
t1=factn2(2*lnpsd+1)/omega
do ias=1,natmtot
  is=idxis(ias)
  do l=0,lmaxo
    t2=t1/(factn2(2*l+1)*rmtl(l,is))
    lma=l**2+1; lmb=lma+2*l
    zlm(lma:lmb)=t2*qlm(lma:lmb,ias)
  end do
! add the pseudocharge and real interstitial densities in G-space
  do ig=1,ngp
    jg=igf(ig)
    if (gpc(ig) > epslat) then
      t2=gpc(ig)*rmt(is)
      t3=1.d0/t2**lnpsd
      z1=t3*fourpi*y00*zlm(1)
      do l=1,lmaxo
        lma=l**2+1; lmb=lma+2*l
        t3=t3*t2
        z1=z1+t3*sum(zlm(lma:lmb)*ylmgp(lma:lmb,ig))
      end do
      z2=jlgprmt(lnpsd,ig,is)*conjg(sfacgp(ig,ias))
      cvclir(jg)=cvclir(jg)+z1*z2
    else
      t2=fourpi*y00/factn2(2*lnpsd+1)
      cvclir(jg)=cvclir(jg)+t2*zlm(1)
    end if
  end do
end do
! solve Poisson's equation in G+p-space for the pseudocharge
do ig=1,ngp
  jg=igf(ig)
  cvclir(jg)=gclgp(ig)*cvclir(jg)
end do
! match potentials at muffin-tin boundary by adding homogeneous solution
do ias=1,natmtot
  is=idxis(ias)
  iro=nri(is)+1
! find the spherical harmonic expansion of the interstitial potential at the
! muffin-tin radius
  zlm(1:lmmaxo)=0.d0
  do ig=1,ngp
    z1=cvclir(igf(ig))*sfacgp(ig,ias)
    zlm(1)=zlm(1)+jlgprmt(0,ig,is)*fourpi*y00*z1
    do l=1,lmaxo
      lma=l**2+1; lmb=lma+2*l
      z2=jlgprmt(l,ig,is)*z1
      zlm(lma:lmb)=zlm(lma:lmb)+z2*conjg(ylmgp(lma:lmb,ig))
    end do
  end do
! add the homogenous solution
  i=np(is)-lmmaxo
  do l=0,lmaxi
    t1=1.d0/rmtl(l,is)
    do lm=l**2+1,(l+1)**2
      z1=t1*(zlm(lm)-cvclmt(i+lm,ias))
      i1=lmmaxi*(nri(is)-1)+lm
      cvclmt(lm:i1:lmmaxi,ias)=cvclmt(lm:i1:lmmaxi,ias)+z1*rl(1:nri(is),l,is)
      i0=i1+lmmaxi
      i1=lmmaxo*(nr(is)-iro)+i0
      cvclmt(i0:i1:lmmaxo,ias)=cvclmt(i0:i1:lmmaxo,ias)+z1*rl(iro:nr(is),l,is)
    end do
  end do
  do l=lmaxi+1,lmaxo
    t1=1.d0/rmtl(l,is)
    do lm=l**2+1,(l+1)**2
      z1=t1*(zlm(lm)-cvclmt(i+lm,ias))
      i0=lmmaxi*nri(is)+lm
      i1=lmmaxo*(nr(is)-iro)+i0
      cvclmt(i0:i1:lmmaxo,ias)=cvclmt(i0:i1:lmmaxo,ias)+z1*rl(iro:nr(is),l,is)
    end do
  end do
end do
! Fourier transform interstitial potential to real-space
call cfftifc(3,ngdg,1,cvclir)
end subroutine

