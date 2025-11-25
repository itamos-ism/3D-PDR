subroutine calc_columndens
use healpix_types
use maincode_module
use global_module
use m_Ray_box
!calculation of column density
!T.Bisbas

#ifdef RAYTHEIA_MO
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(GUIDED) PRIVATE(p,adaptive_step,j,i,k) &
!$OMP PRIVATE(thfpix,phfpix,ray,box1,epray,projected,plength)
#endif
#else
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(GUIDED) PRIVATE(p,adaptive_step,j,i,k)
#endif
#endif
do p=1,pdr_ptot
#ifdef THERMALBALANCE
  if (pdr(p)%fullyconverged) cycle
#endif
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
  if (referee.eq.0) then 
      allocate(pdr(p)%column_NH2(0:nrays-1))
      allocate(pdr(p)%column_NHD(0:nrays-1))
      allocate(pdr(p)%column_NCO(0:nrays-1))
      allocate(pdr(p)%column_NC(0:nrays-1))
      allocate(pdr(p)%column_NS(0:nrays-1))
  endif
  pdr(p)%column_NH2 = 0.0D0
  pdr(p)%column_NHD = 0.0D0
  pdr(p)%column_NCO = 0.0D0
  pdr(p)%column_NC = 0.0D0
  pdr(p)%column_NS = 0.0D0
  do j=0,nrays-1
#ifdef RAYTHEIA_MO
   if (epray(j).eq.0) cycle
   do i=1,epray(j)
      adaptive_step = plength(j,i)
      do k=1,nspec
         if (k.eq.NH2) pdr(p)%column_NH2(j) = pdr(p)%column_NH2(j) + adaptive_step*PC*&
            & pdr(int(projected(j,i)))%rho*pdr(int(projected(j,i)))%abundance(k)

         if (k.eq.NHD) pdr(p)%column_NHD(j) = pdr(p)%column_NHD(j) + adaptive_step*PC*&
            & pdr(int(projected(j,i)))%rho*pdr(int(projected(j,i)))%abundance(k)

         if (k.eq.NCO) pdr(p)%column_NCO(j) = pdr(p)%column_NCO(j) + adaptive_step*PC*&
            & pdr(int(projected(j,i)))%rho*pdr(int(projected(j,i)))%abundance(k)

         if (k.eq.NC) pdr(p)%column_NC(j) = pdr(p)%column_NC(j) + adaptive_step*PC*&
            & pdr(int(projected(j,i)))%rho*pdr(int(projected(j,i)))%abundance(k)

         if (k.eq.NS) pdr(p)%column_NS(j) = pdr(p)%column_NS(j) + adaptive_step*PC*&
            & pdr(int(projected(j,i)))%rho*pdr(int(projected(j,i)))%abundance(k)
      enddo ! End of loop over species (n)
   enddo
#else
     if (pdr(p)%epray(j).eq.0) cycle
     do i=1,pdr(p)%epray(j)!nb_projected_points
#ifdef RAYTHEIA
       adaptive_step = pdr(p)%length(j,i)
       do k=1,nspec
         if (k.eq.NH2) pdr(p)%column_NH2(j) = pdr(p)%column_NH2(j) + adaptive_step*PC*&
            & pdr(int(projected(j,i)))%rho*pdr(int(projected(j,i)))%abundance(k)

         if (k.eq.NHD) pdr(p)%column_NHD(j) = pdr(p)%column_NHD(j) + adaptive_step*PC*&
            & pdr(int(projected(j,i)))%rho*pdr(int(projected(j,i)))%abundance(k)

         if (k.eq.NCO) pdr(p)%column_NCO(j) = pdr(p)%column_NCO(j) + adaptive_step*PC*&
            & pdr(int(projected(j,i)))%rho*pdr(int(projected(j,i)))%abundance(k)

         if (k.eq.NC) pdr(p)%column_NC(j) = pdr(p)%column_NC(j) + adaptive_step*PC*&
            & pdr(int(projected(j,i)))%rho*pdr(int(projected(j,i)))%abundance(k)

         if (k.eq.NS) pdr(p)%column_NS(j) = pdr(p)%column_NS(j) + adaptive_step*PC*&
            & pdr(int(projected(j,i)))%rho*pdr(int(projected(j,i)))%abundance(k)
       enddo ! End of loop over species (n)
#else
       adaptive_step = sqrt((pdr(p)%epoint(1,j,i-1)-pdr(p)%epoint(1,j,i))**2+&
                  &(pdr(p)%epoint(2,j,i-1)-pdr(p)%epoint(2,j,i))**2 + &
                  &(pdr(p)%epoint(3,j,i-1)-pdr(p)%epoint(3,j,i))**2)
       do k=1,nspec
         if (k.eq.NH2) pdr(p)%column_NH2(j) = pdr(p)%column_NH2(j) + adaptive_step*PC*&
             & (pdr(int(pdr(p)%projected(j,i-1)))%rho*pdr(int(pdr(p)%projected(j,i-1)))%abundance(k) +&
             & pdr(int(pdr(p)%projected(j,i)))%rho*pdr(int(pdr(p)%projected(j,i)))%abundance(k))/2.

         if (k.eq.NHD) pdr(p)%column_NHD(j) = pdr(p)%column_NHD(j) + adaptive_step*PC*&
             & (pdr(int(pdr(p)%projected(j,i-1)))%rho*pdr(int(pdr(p)%projected(j,i-1)))%abundance(k) +&
             & pdr(int(pdr(p)%projected(j,i)))%rho*pdr(int(pdr(p)%projected(j,i)))%abundance(k))/2.

         if (k.eq.NCO) pdr(p)%column_NCO(j) = pdr(p)%column_NCO(j) + adaptive_step*PC*&
             & (pdr(int(pdr(p)%projected(j,i-1)))%rho*pdr(int(pdr(p)%projected(j,i-1)))%abundance(k) +&
             & pdr(int(pdr(p)%projected(j,i)))%rho*pdr(int(pdr(p)%projected(j,i)))%abundance(k))/2.

         if (k.eq.NC) pdr(p)%column_NC(j) = pdr(p)%column_NC(j) + adaptive_step*PC*&
             & (pdr(int(pdr(p)%projected(j,i-1)))%rho*pdr(int(pdr(p)%projected(j,i-1)))%abundance(k) +&
             & pdr(int(pdr(p)%projected(j,i)))%rho*pdr(int(pdr(p)%projected(j,i)))%abundance(k))/2.

         if (k.eq.NS) pdr(p)%column_NS(j) = pdr(p)%column_NS(j) + adaptive_step*PC*&
             & (pdr(int(pdr(p)%projected(j,i-1)))%rho*pdr(int(pdr(p)%projected(j,i-1)))%abundance(k) +&
             & pdr(int(pdr(p)%projected(j,i)))%rho*pdr(int(pdr(p)%projected(j,i)))%abundance(k))/2.
       enddo ! End of loop over species (n)
#endif
     enddo ! End of loop over evaluation points along ray (i)
#endif
  enddo ! End of j loop over rays (j)
enddo ! End of ii loop over pdr (ii)
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
return
end subroutine calc_columndens

