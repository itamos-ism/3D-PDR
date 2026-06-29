subroutine opticaldepth
use RTtool_module

!Calculate optical depth
do p=1,itot-1
  step = abs(pdr(p+1)%x-pdr(p)%x)*PC
  do j=1,coo
    do il=1,coolant(j)%cnlev
      do jl=1,coolant(j)%cnlev
        if (jl.ge.il) exit
        if (coolant(j)%A_COEFFS(il,jl).eq.0.0D0) cycle
        sigma = (coolant(j)%FREQUENCIES(il,jl)/C) * SQRT(2*KB*pdr(p)%Tgas/MH/coolant(j)%molweight + vturb**2)
        phi = 1. / sigma 
        frac = 0.5 * ((pdr(p)%coolant(j)%pop(jl)+pdr(p+1)%coolant(j)%pop(jl)) * coolant(j)%WEIGHTS(il)/coolant(j)%WEIGHTS(jl) - &
                (pdr(p)%coolant(j)%pop(il)+pdr(p+1)%coolant(j)%pop(il)))
        if (frac.lt.0.0D0) frac = 0.0D0   ! clip masing lines (no negative opacity)
        pdr(p+1)%coolant(j)%tau(il,jl) = pdr(p)%coolant(j)%tau(il,jl) + phi*(coolant(j)%A_COEFFS(il,jl)*C**2/8./pi/&
                coolant(j)%frequencies(il,jl)**2)*frac*step
      enddo
    enddo
  enddo
enddo

return
end subroutine
