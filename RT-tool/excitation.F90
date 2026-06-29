subroutine excitation
use RTtool_module

!Calculate excitation temperature & black body emission
do p=1,itot
   do j=1,coo
     do il=1,coolant(j)%cnlev
       do jl=1,coolant(j)%cnlev
         if (jl.ge.il) exit
         if (coolant(j)%A_COEFFS(il,jl).eq.0.0D0) cycle
         if (pdr(p)%coolant(j)%pop(jl).eq.0) then
            pdr(p)%coolant(j)%Tex(il,jl) = 0.
            pdr(p)%coolant(j)%Bnu(il,jl) = 0.
         else if (pdr(p)%coolant(j)%pop(jl)*coolant(j)%WEIGHTS(il) .le. &
                 &pdr(p)%coolant(j)%pop(il)*coolant(j)%WEIGHTS(jl)) then
            ! Non-inverted (or equal) populations: the log argument would be
            ! <= 1, giving a non-positive or infinite Tex. Treat as no net line.
            pdr(p)%coolant(j)%Tex(il,jl) = 0.
            pdr(p)%coolant(j)%Bnu(il,jl) = 0.
         else
            pdr(p)%coolant(j)%Tex(il,jl) = (HP*coolant(j)%FREQUENCIES(il,jl)/KB) / log(coolant(j)%WEIGHTS(il)*&
                    pdr(p)%coolant(j)%pop(jl)/pdr(p)%coolant(j)%pop(il)/coolant(j)%WEIGHTS(jl))
            ! Guard the Tex singularity as the population ratio -> 1+. Tex cannot
            ! physically exceed the local gas temperature for these collisionally
            ! excited lines, so clamp it to Tgas.
            if (pdr(p)%coolant(j)%Tex(il,jl) .gt. pdr(p)%Tgas) &
                pdr(p)%coolant(j)%Tex(il,jl) = pdr(p)%Tgas
            pdr(p)%coolant(j)%Bnu(il,jl) = (2.*HP*coolant(j)%FREQUENCIES(il,jl)**3/C**2) / (dexp(HP*&
                    coolant(j)%FREQUENCIES(il,jl)/KB/pdr(p)%coolant(j)%Tex(il,jl))-1.0)
            Bckgr = (2.*HP*coolant(j)%FREQUENCIES(il,jl)**3/C**2) / (dexp(HP*coolant(j)%FREQUENCIES(il,jl)/KB/Tcmb)-1.0)

            pdr(p)%coolant(j)%Bnu(il,jl) = pdr(p)%coolant(j)%Bnu(il,jl) - Bckgr
          endif
        enddo
      enddo
    enddo
enddo

return
end subroutine
