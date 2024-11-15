subroutine bfuinit
use modmain
use modulr
use modrandom
!
implicit none
!
real(8) cb,t1
integer ir, idm, ias,ia,is
real(8) rfft(nqpt)
!
real(8),allocatable ::  bfcr(:,:),bfcmtr(:,:,:)
!
allocate(bfcr(3,nqpt),bfcmtr(natmtot,3,nqpt))
bfcr(:,:)=0.d0
bfcmtr(:,:,:)=0.d0
if (trdbfcr) then
! if read-in the R-dependent magnetic field from external file.
  call readpotvecr('BFCR.OUT',bfcr(1:3,1:nqpt))
  write(*,'("Info(bfuinit) : external magnetic field read in from file ",A)')'BFCR.OUT'
else    
  do ir=1,nqpt
    do idm=1,3
      t1=0.d0
! Apply a random magnetic field if the magnitude is bigger than 1( 1.d-7 ~= 1.715) Gauss
      if (rndbfcu>=1.d-7) t1=rndbfcu*(randomu()-0.5d0)
      bfcr(idm,ir)=bfieldc(idm)+t1
! Not apply random magnetic field to the local muffin-tin 
      do ias=1,natmtot
        is=idxis(ias)
        ia=idxia(ias)
        bfcmtr(ias,idm,ir)=bfcmt(idm,ia,is)
      end do
    end do
  end do
! Output the R-dependent magnetic field to the file.
  call writepotvecr('BFCR.OUT', bfcr(1:3,1:nqpt))
end if
! Prefactor
cb=gfacte/(4.d0*solsc)
! Fourier transform R-dependent magnetic filed to Q-mesh and multiplied by the prefactor
rfft(1:nqpt)=cb*bfcr(3,1:nqpt)/omega
call rzfftifc(3,ngridq,-1,rfft,bfcq(ndmag,1:nfqrz))
do ias=1,natmtot
  rfft(1:nqpt)=cb*bfcmtr(ias,3,1:nqpt)/omega
  call rzfftifc(3,ngridq,-1,rfft(1:nqpt),bfcmtq(ias,ndmag,1:nfqrz))
end do
if (ncmag) then
  do idm=1,2
    rfft(1:nqpt)=cb*bfcr(idm,1:nqpt)/omega
    call rzfftifc(3,ngridq,-1,rfft,bfcq(idm,1:nfqrz))
    do ias=1,natmtot
      rfft(1:nqpt)=cb*bfcmtr(ias,idm,1:nqpt)/omega
      call rzfftifc(3,ngridq,-1,rfft(1:nqpt),bfcmtq(ias,idm,1:nfqrz))
    end do
  end do
end if
bfcq(1:ndmag,1)=dble(bfcq(1:ndmag,1))
bfcmtq(1:natmtot,1:ndmag,1)=dble(bfcmtq(1:natmtot,1:ndmag,1))
deallocate(bfcr,bfcmtr)
end subroutine

