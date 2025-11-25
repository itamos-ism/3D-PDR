! Calculation of the UVfield
SUBROUTINE CALC_UVFIELD

use healpix_types
use maincode_module
use global_module
#ifdef RAYTHEIA_MO
use m_Mesh
use m_Ray_box
#else
#ifdef RAYTHEIA
use m_Mesh
use m_Ray_box
#endif
#endif

real(kind=dp) :: UVdot
real(kind=dp) :: UVdotsum
integer(kind=I4B) :: UVinter

do p=1,pdr_ptot
      pdr(p)%UVfield = 0.0D0
enddo

#ifdef RAYTHEIA_MO
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,adaptive_step,j,i) &
!$OMP PRIVATE(thfpix,phfpix,ray,box1,epray,projected,plength,UVdot,UVdotsum,UVinter)
#endif
#else
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,adaptive_step,j,i)
#endif
#endif
DO p=1,pdr_ptot
  pdr(p)%rad_surface = 0.0D0
  pdr(p)%AV = 0.0D0
#ifdef RAYTHEIA_MO
   box1%min = corner_min
   box1%max = corner_max
   epray = 0
   projected = p
   plength = 0.D0
   do j=0,nrays-1
      call pix2ang_nest(nside, j, thfpix, phfpix)
      ray%origin = [pdr(p)%x, pdr(p)%y, pdr(p)%z]
      ray%angle = [thfpix, phfpix]
      call raytheia(ray, box1, levels-1, j, epray, projected, plength)
   enddo
#else
  pdr(p)%projected(:,0) = p
#endif
  DO j=0,nrays-1
#ifdef RAYTHEIA_MO
     IF (epray(j).GT.0) THEN
         do i=1,epray(j)
            adaptive_step = plength(j,i)
            pdr(p)%AV(j) = pdr(p)%AV(j) + &
               &pdr(INT(projected(j,i)))%rho*adaptive_step*pc      
         enddo
     ENDIF
#else
#ifdef RAYTHEIA
     IF (pdr(p)%epray(j).GT.0) THEN
        DO i=1,pdr(p)%epray(j)
           adaptive_step = pdr(p)%length(j,i)
           pdr(p)%AV(j) = pdr(p)%AV(j) + &
              &pdr(INT(pdr(p)%projected(j,i)))%rho*adaptive_step*pc
        ENDDO
     ENDIF
#else     
     IF (pdr(p)%epray(j).GT.0) THEN
        DO i=1,pdr(p)%epray(j)
           adaptive_step = SQRT((pdr(p)%epoint(1,j,i-1)-pdr(p)%epoint(1,j,i))**2 + &
                              & (pdr(p)%epoint(2,j,i-1)-pdr(p)%epoint(2,j,i))**2 + &
                              & (pdr(p)%epoint(3,j,i-1)-pdr(p)%epoint(3,j,i))**2)
           pdr(p)%AV(j) = pdr(p)%AV(j) + &
              &((pdr(INT(pdr(p)%projected(j,i-1)))%rho + &
              &pdr(INT(pdr(p)%projected(j,i)))%rho)/2.)*adaptive_step*pc
        ENDDO
     ENDIF
#endif
#endif
     pdr(p)%AV(j) = pdr(p)%AV(j)*AV_fac
     if (pdr(p)%AV(j).gt.100.0) pdr(p)%AV(j)=100.0D0 !limiting AV to avoid numerical errors
  ENDDO
#ifdef ONEDIMENSIONAL
   pdr(p)%rad_surface = 0.0D0
   pdr(p)%rad_surface(6) = Gext(1)!-DOT_PRODUCT(Gext(:),vectors(:,j))
#else
   pdr(p)%rad_surface = Gext(1)/real(nrays,kind=dp)
#endif

   DO j=0,nrays-1
      pdr(p)%UVfield = pdr(p)%UVfield + pdr(p)%rad_surface(j)*EXP(-pdr(p)%AV(j)*UV_fac)
      IF (pdr(p)%UVfield.LT.1.0D-50) pdr(p)%UVfield = 0.0D0
   ENDDO ! End of loop over rays
ENDDO ! End of loop over pdr cells
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef THERMALBALANCE
#ifdef GUESS_TEMP
DO p=1,pdr_ptot
   Tguess = 10.0D0*(1.0D0+(2.*pdr(p)%UVfield)**(1.0D0/3.0D0))
   pdr(p)%nTgas = Tguess
   pdr(p)%Tgas = Tguess

   pdr(p)%Tlow = Tguess/2.0D0
   pdr(p)%Thigh = Tguess*1.5D0

   IF (pdr(p)%Tlow.LT.Tmin)  pdr(p)%Tlow  = Tmin
   IF (pdr(p)%Thigh.GT.Tmax) pdr(p)%Thigh = Tmax
ENDDO
#endif
#endif
RETURN

END SUBROUTINE CALC_UVFIELD
