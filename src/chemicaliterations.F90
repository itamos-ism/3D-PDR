subroutine chemicaliterations(iw,ichem)
use healpix_types
use maincode_module
use global_module
use maincode_local
integer,intent(in)::iw,ichem
integer,parameter::npart=1

DO II=1,ichem
  write(6,'(" Chemical iteration = ",I2)') ii

#ifdef OPENMP
!$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(SHARED) PRIVATE (p,rate,NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
#endif
   do p=1,pdr_ptot
#ifdef THERMALBALANCE
       if (iw.eq.2.and.pdr(p)%fullyconverged) cycle
#endif
       CALL CALCULATE_REACTION_RATES(pdr(p)%nTgas,pdr(p)%Tdust,nrays,pdr(p)%rad_surface(0:nrays-1),pdr(p)%AV(0:nrays-1),&
             &pdr(p)%column_NH2(0:nrays-1),pdr(p)%column_NHD(0:nrays-1),pdr(p)%column_NCO(0:nrays-1),&
             &pdr(p)%column_NC(0:nrays-1),pdr(p)%column_NS(0:nrays-1),&
             &nreac, reactant, product, alpha, beta, gamma, rate, rtmin, rtmax, duplicate, nspec,&
             &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI,pdr(p)%abundance(NELECT)*pdr(p)%rho,pdr(p)%rho,pdr(p)%zetalocal)
       call calculate_abundances(pdr(p)%abundance,rate,pdr(p)%rho,pdr(p)%nTgas,npart,nspec,nreac)
#ifdef CHEMANALYSIS
       temp_rate(:,p) = rate
#endif
   enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

   call calc_columndens
ENDDO 

#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,i,temp_Z_function) 
#endif
do p=1,pdr_ptot
  do i=1,coo
    call  calculate_partition_function(temp_Z_function,coolant(i)%cnlev,&
            coolant(i)%energies,coolant(i)%weights,pdr(p)%nTgas)
    call calculate_lte_populations(coolant(i)%cnlev,pdr(p)%coolant(i)%pop,coolant(i)%energies,&
            coolant(i)%weights,temp_Z_function,pdr(p)%abundance(coolant(i)%cspec)*pdr(p)%rho,&
            pdr(p)%nTgas)
  enddo
enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif


return
end subroutine
