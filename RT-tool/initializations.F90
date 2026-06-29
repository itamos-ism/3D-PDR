subroutine initializations
use RTtool_module

allocate(Nabn(1:nspec))
pdrfile = trim(adjustl(filein))//".pdr.fin"
close(1);open(unit=1,file=pdrfile,status='old')
itot=0
do 
  read(1,*,end=100) dummy
  itot=itot+1
enddo
100 continue
close(1)

allocate(pdr(0:itot))
allocate(pdr(0)%coolant(1:coo))
do p=1,itot
  allocate(pdr(p)%Nabn(1:nspec))
  allocate(pdr(p)%coolant(1:coo))
  do j=1,coo
    allocate(pdr(p)%coolant(j)%tau(1:coolant(j)%cnlev,1:coolant(j)%cnlev))
    allocate(pdr(p)%coolant(j)%Tr(1:coolant(j)%cnlev,1:coolant(j)%cnlev))
    allocate(pdr(p)%coolant(j)%Tex(1:coolant(j)%cnlev,1:coolant(j)%cnlev))
    allocate(pdr(p)%coolant(j)%Bnu(1:coolant(j)%cnlev,1:coolant(j)%cnlev))
    pdr(p)%coolant(j)%tau=0.0D0
    pdr(p)%coolant(j)%Tr=0.0D0
    pdr(p)%coolant(j)%Tex=0.0D0
    pdr(p)%coolant(j)%Bnu=0.0D0
  enddo
enddo
do j=1,coo
  allocate(pdr(0)%coolant(j)%tau(1:coolant(j)%cnlev,1:coolant(j)%cnlev))
  allocate(coolant(j)%Tr(1:coolant(j)%cnlev,1:coolant(j)%cnlev))
  allocate(coolant(j)%units(1:coolant(j)%cnlev,1:coolant(j)%cnlev))
  coolant(j)%Tr=0.0D0
  coolant(j)%N=0.0D0
  pdr(0)%coolant(j)%tau=0.0D0
  do il=1,coolant(j)%cnlev
    do jl=1,coolant(j)%cnlev
      if (jl.ge.il) exit
      if (coolant(j)%A_COEFFS(il,jl).eq.0.0D0) cycle
      coolant(j)%units(il,jl) = C**2/2./KB/coolant(j)%frequencies(il,jl)**2
    enddo
  enddo
enddo

return
end subroutine
