subroutine initialization
use maincode_module

do p=1,pdr_ptot
  allocate(pdr(p)%abundance(1:nspec))
  pdr(p)%abundance = init_abundance
  pdr(p)%levelconverged = .false.
  pdr(p)%previouschange = 'N'
  pdr(p)%totalcooling = 0.0D0
  allocate(pdr(p)%cooling(coo), pdr(p)%heating(12))
  pdr(p)%cooling = 0.0D0
  pdr(p)%heating = 0.0D0
#ifdef THERMALBALANCE
  ! The first column calculation precedes the density-based freeze below.
  pdr(p)%fullyconverged = .false.
  pdr(p)%doleveltmin = .false.
  pdr(p)%Fmean = 0.0D0
  pdr(p)%Fratio = 0.0D0
  pdr(p)%dobinarychop = .false.
#ifdef ILLINOIS
  pdr(p)%ill_last = "N"
  pdr(p)%Flow = 0.0D0
  pdr(p)%Fhigh = 0.0D0
#endif
#ifdef REBRACKET
  pdr(p)%nrebracket = 0
#endif
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
    pdr(p)%coolant(i)%pop = 0.0D0
    pdr(p)%coolant(i)%solution = 0.0D0
    pdr(p)%coolant(i)%line = 0.0D0
    pdr(p)%coolant(i)%relativechange = 0.0D0
    pdr(p)%coolant(i)%isconverged = .false.
#ifdef NGACCEL
    allocate(pdr(p)%coolant(i)%pophist(coolant(i)%cnlev,3))
    pdr(p)%coolant(i)%pophist = 0.0D0
#endif
#ifdef NGRELAX
    pdr(p)%coolant(i)%noscil = 0
#endif
    pdr(p)%coolant(i)%coolprev = 0.0D0
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
