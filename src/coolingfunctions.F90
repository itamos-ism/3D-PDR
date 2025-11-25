subroutine coolingfunctions
use healpix_types
use maincode_module
use global_module
use m_Ray_box

#ifdef RAYTHEIA_MO
#ifdef OPENMP
!$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(SHARED) PRIVATE(cpop,temp_line,p,j,i,k,ilevel) &
!$OMP PRIVATE(temp_C_COEFFS, temp_transition, temp_solution) &
!$OMP PRIVATE(thfpix,phfpix,ray,box1,epray,projected,plength)
#endif
#else
#ifdef OPENMP
!$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(SHARED) PRIVATE(cpop,temp_line,p,j,i,k,ilevel) & 
!$OMP PRIVATE(temp_C_COEFFS, temp_transition, temp_solution)
#endif
#endif
do p=1,pdr_ptot
!----------------------------------------------
#ifdef THERMALBALANCE
      if (pdr(p)%levelconverged.or.pdr(p)%fullyconverged) cycle
#else
      if (pdr(p)%levelconverged) cycle
#endif
!----------------------------------------------
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

      ! Specify the evaluation points along each ray from the current pdrpoint
      allocate(cpop(1:coo))
      do k = 1, coo
          allocate(cpop(k)%evalpop(0:nrays-1,0:maxpoints,1:coolant(k)%cnlev))
          cpop(k)%evalpop=0.0D0
          do j = 0, nrays - 1
#ifdef RAYTHEIA_MO
            do i = 0, epray(j)
                do ilevel = 1, coolant(k)%cnlev
                   cpop(k)%evalpop(j,i,ilevel) = pdr(int(projected(j,i)))%coolant(k)%pop(ilevel)
                enddo
             enddo
          enddo
#else
            do i = 0, pdr(p)%epray(j)
                do ilevel = 1, coolant(k)%cnlev
                   cpop(k)%evalpop(j,i,ilevel) = pdr(int(pdr(p)%projected(j, i)))%coolant(k)%pop(ilevel)
                enddo
             enddo
          enddo
#endif
    
          ! Use the LVG (escape probability) method to determine the transition matrices and solve for the level populations
          allocate(temp_C_COEFFS(1:coolant(k)%cnlev,1:coolant(k)%cnlev))
          allocate(temp_line(1:coolant(k)%cnlev,1:coolant(k)%cnlev))
          allocate(temp_transition(1:coolant(k)%cnlev,1:coolant(k)%cnlev))
          allocate(temp_solution(1:coolant(k)%cnlev))
          CALL FIND_CCOEFF(coolant(k)%cntemp,coolant(k)%cnlev,pdr(p)%nTgas,coolant(k)%temperatures,&
             & coolant(k)%H_COL,coolant(k)%HP_COL,coolant(k)%EL_COL,coolant(k)%HE_COL,coolant(k)%H2_COL,& 
             & coolant(k)%PH2_COL,coolant(k)%OH2_COL,temp_C_COEFFS,pdr(p)%abundance(NH)*pdr(p)%rho,&
             & pdr(p)%abundance(NPROTON)*pdr(p)%rho, pdr(p)%abundance(NELECT)*pdr(p)%rho, &
             & pdr(p)%abundance(NHE)*pdr(p)%rho,pdr(p)%abundance(NH2)*pdr(p)%rho)
#ifdef RAYTHEIA_MO
          call escape_probability(temp_transition, pdr(p)%Tdust, nrays, coolant(k)%cnlev, &
                  &coolant(k)%A_COEFFS, coolant(k)%B_COEFFS, temp_C_COEFFS, &
                  &coolant(k)%frequencies, cpop(k)%evalpop, maxpoints, &
                  &pdr(p)%nTgas, v_turb, epray, pdr(p)%coolant(k)%pop, &
                  &coolant(k)%weights,pdr(p)%cooling(k),temp_line,&
                  &pdr(p)%rho,metallicity,coolant(k)%molweight, plength)
#else
#ifdef RAYTHEIA
          call escape_probability(temp_transition, pdr(p)%Tdust, nrays, coolant(k)%cnlev, &
                  &coolant(k)%A_COEFFS, coolant(k)%B_COEFFS, temp_C_COEFFS, &
                  &coolant(k)%frequencies, cpop(k)%evalpop, maxpoints, &
                  &pdr(p)%nTgas, v_turb, pdr(p)%epray, pdr(p)%coolant(k)%pop, &
                  &coolant(k)%weights,pdr(p)%cooling(k),temp_line,&
                  &pdr(p)%rho,metallicity,coolant(k)%molweight,pdr(p)%length)
#else
          call escape_probability(temp_transition, pdr(p)%Tdust, nrays, coolant(k)%cnlev, &
                 &coolant(k)%A_COEFFS, coolant(k)%B_COEFFS, temp_C_COEFFS, &
                 &coolant(k)%frequencies, cpop(k)%evalpop, maxpoints, &
                 &pdr(p)%nTgas, v_turb, pdr(p)%epray, pdr(p)%coolant(k)%pop, &
                 &pdr(p)%epoint,coolant(k)%weights,pdr(p)%cooling(k),temp_line,&
                 &pdr(p)%rho,metallicity,coolant(k)%molweight)
#endif
#endif
          pdr(p)%coolant(k)%line=temp_line
          call solvlevpop(coolant(k)%cnlev,temp_transition,pdr(p)%abundance(coolant(k)%cspec)*pdr(p)%rho,&
                  temp_solution)
          pdr(p)%coolant(k)%solution=temp_solution
          deallocate(temp_C_COEFFS)
          deallocate(temp_line)
          deallocate(temp_transition)
          deallocate(temp_solution)
          deallocate(cpop(k)%evalpop)
      enddo
      deallocate(cpop)

#ifdef FORCECONVERGENCE
      ! If the level populations are oscillating and not converging, try to suppress the oscillations
      ! by taking the average of the level populations from this iteration and the previous iteration
      if (all([(int(coolant(k)%percentage)==100, k=2,4)])) then !special treatment of coolants other than C+,C,O
           if (levpop_iteration.ge.120) then
                   pdr(p)%coolant(1)%solution=pdr(p)%coolant(1)%pop
           else if (levpop_iteration.ge.75) then
                   pdr(p)%coolant(1)%solution=0.5*(pdr(p)%coolant(1)%solution + pdr(p)%coolant(1)%pop)
           endif
      endif
      do k=2,coo
           if (levpop_iteration.ge.120) then
               pdr(p)%coolant(k)%solution=pdr(p)%coolant(k)%pop
           else if (levpop_iteration.ge.75) then
               pdr(p)%coolant(k)%solution=0.5*(pdr(p)%coolant(k)%solution + pdr(p)%coolant(k)%pop)
           endif
      enddo
#endif
            !-------------------------------------------------------------------

      pdr(p)%totalcooling = sum(pdr(p)%cooling(:))

enddo !particles
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

return
end subroutine
