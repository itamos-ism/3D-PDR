subroutine readpdrfiles
use RTtool_module

close(1);open(unit=1,file='chemfiles/species_'//trim(adjustl(chemsuf)),status='old')
do i=1,nspec
  read(1,*) index,species(i)
  do j=1,coo
    if (coolant(j)%cname.eq.species(i)) coolant(j)%cspec=i
  enddo
enddo

pdrfile = trim(adjustl(filein))//".pdr.fin"

close(1);open(unit=1,file=pdrfile,status='old')

do p=1,itot
  allocate(pdr(p)%abun(1:nspec))
  read(1,*) id,pdr(p)%x,pdr(p)%Av,pdr(p)%Tgas,pdr(p)%Tdust,etype,pdr(p)%rho,uv,pdr(p)%abun
enddo
do p=1,itot
  do j=1,coo
    allocate(pdr(p)%coolant(j)%pop(1:coolant(j)%cnlev))
  enddo
enddo

return
end subroutine

