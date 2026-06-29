subroutine writeoutputs
use RTtool_module


write(6,*) '' 
write(6,*) '<n(H2)> = ',avg_nH2
write(6,*) '<Tgas> = ',avg_Tgas
write(6,*) 'Tcmb = ',Tcmb

write(6,*) ''
do j=1,coo
  write(6,'("N(",A,")=",ES15.7," cm-2")') trim(adjustl(coolant(j)%cname)),coolant(j)%N
  write(6,*) '--------------------'
  write(6,'(A10,A15,A15,A15,A15)') 'Line', 'Freq.(GHz)', 'Tex (K)', 'tau', 'Tr (K)'
  if (coolant(j)%cname=='C+') coolev=1
  if (coolant(j)%cname=='C') coolev=2
  if (coolant(j)%cname=='O') coolev=2
  if (coolant(j)%cnlev>10) coolev=10
  do il=1,coolant(j)%cnlev
    do jl=1,coolant(j)%cnlev
      if (jl.ge.il) exit
      if (coolant(j)%A_COEFFS(il,jl).eq.0.0D0) cycle
        write(6,'(A,"(",I2,"-",I2,") :",4ES15.7)') trim(coolant(j)%cname), il-1, jl-1, &
          coolant(j)%frequencies(il,jl)/1d9, &
          pdr(itot)%coolant(j)%Tex(il,jl),pdr(itot-1)%coolant(j)%tau(il,jl), coolant(j)%Tr(il,jl)*coolant(j)%units(il,jl)
    enddo
  enddo
  write(6,*) ''
enddo
    
do j=1,coo
  close(1);open(unit=1,file=adjustl(trim(directory))//'/Tr_'//adjustl(trim(prefix))//'.'//&
          adjustl(trim(coolant(j)%cname))//'.dat',status='replace')
  close(2);open(unit=2,file=adjustl(trim(directory))//'/tau_'//adjustl(trim(prefix))//'.'//&
          adjustl(trim(coolant(j)%cname))//'.dat',status='replace')
  close(4);open(unit=4,file=adjustl(trim(directory))//'/Tex_'//adjustl(trim(prefix))//'.'//&
          adjustl(trim(coolant(j)%cname))//'.dat',status='replace')
  if (coolant(j)%cname=='C+') coolev=1
  if (coolant(j)%cname=='C') coolev=2
  if (coolant(j)%cname=='O') coolev=2
  if (coolant(j)%cnlev>10) coolev=10

  write(1,*) '# Ntot, NH2, Tgas, Av, Nmol, Tr...'
  write(4,*) '# Tex... (row-aligned with Tr/tau)'
  do p=1,itot-1
    !all outputs contain Ntot, N(H2), Tgas, Av, N(coolant) per depth point
    write(1,'(5ES15.7)',advance='no') pdr(p)%Ntot,pdr(p)%Nabn(iH2),pdr(p)%Tgas,pdr(p)%Av,pdr(p)%coolant(j)%N
    do l=1,coolev
       write(1,'(100ES15.7)',advance='no') pdr(p)%coolant(j)%Tr(l+1,l)*coolant(j)%units(l+1,l)
       write(2,'(100ES15.7)',advance='no') pdr(p)%coolant(j)%tau(l+1,l)
       write(4,'(100ES15.7)',advance='no') pdr(p)%coolant(j)%Tex(l+1,l)
    enddo
    write(1,*);write(2,*);write(4,*)
  enddo
enddo

write(33,*) '#',coolant(:)%cname
do p=1,itot-1
  do j=1,coo
    write(33,'(100ES15.7)',advance='no') pdr(p)%coolant(j)%N
  enddo
  write(33,*)
enddo


return
end subroutine
