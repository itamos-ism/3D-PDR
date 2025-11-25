subroutine readcoolants
use maincode_module

allocate(coolant(1:coo))
do i = 1,coo
  open(unit=44,file=coolfile(i),status='old')
  read(44,'()') 
  read(44,*) coolant(i)%cname
  read(44,'()')
  read(44,*) coolant(i)%molweight
  read(44,'()')
  read(44,*) cur_nlev, cur_ntemp
  close(44)
  coolant(i)%cnlev = cur_nlev
  coolant(i)%cntemp = cur_ntemp
  allocate(coolant(i)%energies(1:cur_nlev))
  allocate(coolant(i)%weights(1:cur_nlev))
  allocate(coolant(i)%A_COEFFS(1:cur_nlev,1:cur_nlev))
  allocate(coolant(i)%B_COEFFS(1:cur_nlev,1:cur_nlev))
  allocate(coolant(i)%frequencies(1:cur_nlev,1:cur_nlev))
  allocate(coolant(i)%temperatures(1:7,1:cur_ntemp))
  allocate(coolant(i)%H_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%HP_COL(1:cur_nlev,1:cur_nlev,1:cur_ntemp))
  allocate(coolant(i)%EL_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%HE_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%H2_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%PH2_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%OH2_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  call readinput(coolfile(i),cur_nlev,cur_ntemp,coolant(i)%energies,coolant(i)%weights,&
          coolant(i)%A_COEFFS,coolant(i)%B_COEFFS,coolant(i)%frequencies,coolant(i)%temperatures,&
          coolant(i)%H_COL,coolant(i)%HP_COL,coolant(i)%EL_COL,coolant(i)%HE_COL,coolant(i)%H2_COL,&
          coolant(i)%PH2_COL,coolant(i)%OH2_COL)
enddo
return
end subroutine
