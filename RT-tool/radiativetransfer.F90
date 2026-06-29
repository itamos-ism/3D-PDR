subroutine radiativetransfer
use RTtool_module

!Solve radiative transfer equation
Ntot = 0; Nabn = 0
Ntgas = 0; Nnh2 = 0

do p=1,itot-1
   step = abs(pdr(p+1)%x-pdr(p)%x)*PC
   Ntot = Ntot + 0.5*(pdr(p)%rho+pdr(p+1)%rho)*step
   Nabn = Nabn + 0.5*(pdr(p)%rho*pdr(p)%abun + pdr(p+1)%rho*pdr(p+1)%abun)*step
   Ntgas = Ntgas + 0.5*(pdr(p)%rho*pdr(p)%Tgas + pdr(p+1)%rho*pdr(p+1)%Tgas)*step
   Nnh2 = Nnh2 + 0.5*(pdr(p)%rho**2*pdr(p)%abun(iH2) + pdr(p+1)%rho**2*pdr(p+1)%abun(iH2))*step
   pdr(p)%Ntot = Ntot
   pdr(p)%Nabn = Nabn

   do j=1,coo
      coolant(j)%N = coolant(j)%N + 0.5*(pdr(p)%rho*pdr(p)%abun(coolant(j)%cspec) + &
                        pdr(p+1)%rho*pdr(p+1)%abun(coolant(j)%cspec))*step
      pdr(p)%coolant(j)%N = coolant(j)%N !stores increment
      do il=1,coolant(j)%cnlev
        do jl=1,coolant(j)%cnlev
          if (jl.ge.il) exit
          if (coolant(j)%A_COEFFS(il,jl).eq.0.0D0) cycle
          dtau=abs(pdr(p)%coolant(j)%tau(il,jl)-pdr(p+1)%coolant(j)%tau(il,jl))
          if (dtau.gt.1d10) then
             coolant(j)%Tr(il,jl) = pdr(p)%coolant(j)%Bnu(il,jl)
          else if (dtau.gt.1d-6) then
             Tr_incr1 = pdr(p)%coolant(j)%Bnu(il,jl) * ((1.-dexp(-dtau))/dtau-dexp(-dtau))
             Tr_incr2 = pdr(p+1)%coolant(j)%Bnu(il,jl) * (1.-(1.-dexp(-dtau))/dtau)
             coolant(j)%Tr(il,jl) = coolant(j)%Tr(il,jl) * dexp(-dtau) + Tr_incr1 + Tr_incr2
          else
             coolant(j)%Tr(il,jl) = coolant(j)%Tr(il,jl) * (1.-dtau) + (pdr(p)%coolant(j)%Bnu(il,jl) + &
                     pdr(p+1)%coolant(j)%Bnu(il,jl))*dtau/2.
          endif

          pdr(p)%coolant(j)%Tr(il,jl) = coolant(j)%Tr(il,jl)
          if (j.eq.coo.and.il.eq.3) write(22,*) il,jl,pdr(p)%x,pdr(p)%coolant(j)%Tr(il,jl),coolant(j)%units(il,jl)
        enddo
      enddo
   enddo
enddo

avg_nH2 = Nnh2 / Ntot
avg_Tgas = Ntgas / Ntot


return
end subroutine
