!Calculate CR attenuation
subroutine calc_CRattenuation

use healpix_types
use maincode_module
use global_module

implicit none
integer :: j
real(kind=dp) :: coltest, crfieldaux

#ifdef CRATTENUATION

DO p=1,pdr_ptot
  pdr(p)%zetalocal = 0d0
#ifndef ONEDIMENSIONAL
  coltest = 0d0
  DO j=0,nrays-1
    coltest = coltest + EXP(-3.5*pdr(p)%AV(j))
  ENDDO
  coltest = (-1.0/3.5)*LOG((1.0/real(nrays,kind=dp))*coltest)
  coltest = coltest/AV_fac
  coltest = max(coltest, 5e17)
#if CRATTENUATION == 1
  call zetafuncPad(coltest, pdr(p)%zetalocal)
#else 
  call zetafuncPoly(coltest, pdr(p)%zetalocal)
#endif
  pdr(p)%zetalocal = pdr(p)%zetalocal
#else
#if CRATTENUATION == 1
  call zetafuncPad(pdr(p)%AV(6)/AV_fac, pdr(p)%zetalocal)
#else
  call zetafuncPoly(pdr(p)%AV(6)/AV_fac, pdr(p)%zetalocal)
#endif
#endif
  pdr(p)%zetalocal = pdr(p)%zetalocal/1.3d-17
enddo

#else

do p=1,pdr_ptot
  pdr(p)%zetalocal = zeta
end do

#endif

return 

contains

subroutine zetafuncPad(ncol, zetalocal)
  use definitions
  implicit none
  real(kind=DP), intent(in) :: ncol
  real(kind=DP), intent(out):: zetalocal
  integer :: ii

  real(kind=DDP), dimension(0:9) :: fl = (/ -3.331056497233d6, &
  & 1.207744586503d6, -1.913914106234d5, 1.731822350618d4, &
  & -9.790557206178d2, 3.543830893824d1, -8.03486945420d-01, &
  & 1.048808593086d-02, -6.188760100997d-05, 3.122820990797d-08/)

  real(kind=DDP), dimension(0:9) :: fh = (/ 1.001098610761d7, &
  & -4.231294690194d6, 7.921914432011d5, -8.623677095423d4, &
  & 6.015889127529d3, -2.789238383353d2, 8.595814402406d0, &
  & -1.698029737474d-1, 1.951179287567d-3, -9.937499546711d-6 /)

  real(kind=DDP) :: craux

  real(kind=DDP) :: lcol

  craux = 0d0
  lcol = LOG10(MAX(ncol, 5E18))

  do ii=0,9
    if ((crfieldchoice .eq. "L") .or. (crfieldchoice .eq. "l")) then
      craux = craux + fl(ii)*lcol**ii
    endif
    if ((crfieldchoice .eq. "H") .or. (crfieldchoice .eq. "h")) then
      craux = craux + fh(ii)*lcol**ii
    endif
  end do

  zetalocal = 10**craux

return 
end subroutine

subroutine zetafuncPoly(ncol, zetalocal)
  use definitions
  use maincode_module, only: crattennorm, crattenslope, crattenn0
  implicit none
  real(kind=DP), intent(in) :: ncol
  real(kind=DP), intent(out):: zetalocal
  real(kind=DP) :: craux, critsurface

  critsurface = 4.182e25 ! This is 98 g/cm^3
  zetalocal = crattennorm * ((1.0 + ncol/crattenn0)**crattenslope)* EXP(-1.0*ncol/critsurface)

end subroutine

end subroutine
