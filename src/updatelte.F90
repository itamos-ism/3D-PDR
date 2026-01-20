subroutine updatelte

use healpix_types
use maincode_module
use global_module
use maincode_local

#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,i,temp_Z_function) 
#endif
do p=1,pdr_ptot
#ifdef THERMALBALANCE
  if (pdr(p)%fullyconverged) cycle
#endif
  do i=2,coo
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
