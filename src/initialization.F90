subroutine initialization
use maincode_module

do p=1,pdr_ptot
  allocate(pdr(p)%abundance(1:nspec))
  pdr(p)%abundance = init_abundance
#ifdef THERMALBALANCE
  pdr(p)%dobinarychop = .false.
#endif
#ifdef RESTART
  pdr(p)%restconverged = .false.
#endif
enddo

write(6,*) 'Initialization'
#ifdef RAYTHEIA_MO
  allocate(epray(0:nrays-1))
  allocate(plength(0:nrays-1,0:maxpoints))
  allocate(projected(0:nrays-1,0:maxpoints))
#endif
do p=1,pdr_ptot
  allocate(pdr(p)%coolant(1:coo))
  do i=1,coo
    allocate(pdr(p)%coolant(i)%pop(coolant(i)%cnlev)) 
    allocate(pdr(p)%coolant(i)%line(coolant(i)%cnlev,coolant(i)%cnlev))
    allocate(pdr(p)%coolant(i)%solution(coolant(i)%cnlev))
    allocate(pdr(p)%coolant(i)%relativechange(coolant(i)%cnlev))
  enddo
#ifndef RAYTHEIA_MO
  allocate(pdr(p)%epray(0:nrays-1))                   
#ifndef RAYTHEIA
  allocate(pdr(p)%epoint(1:3,0:nrays-1,0:maxpoints))  
#else
  allocate(pdr(p)%length(0:nrays-1,0:maxpoints))
#endif
  allocate(pdr(p)%projected(0:nrays-1,0:maxpoints))
#endif   
  allocate(pdr(p)%raytype(0:nrays-1))                 
  allocate(pdr(p)%AV(0:nrays-1))
  allocate(pdr(p)%rad_surface(0:nrays-1))
  allocate(pdr(p)%column_NH2(0:nrays-1))
  allocate(pdr(p)%column_NHD(0:nrays-1))
  allocate(pdr(p)%column_NCO(0:nrays-1))
  allocate(pdr(p)%column_NC(0:nrays-1))
  allocate(pdr(p)%column_NS(0:nrays-1))
enddo

#ifndef GUESS_TEMP
do p=1,pdr_ptot
  pdr(p)%nTgas = Tguess
  pdr(p)%Tgas = Tguess
#ifdef THERMALBALANCE
  pdr(p)%Tlow = Tlow0
  pdr(p)%Thigh = Thigh0
#endif
enddo
#endif

return
end subroutine
