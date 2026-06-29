subroutine readspopfile
use RTtool_module

!reading RTspop.fin header
write(6,*) 'Reading .RTspop.fin'
open(unit=2,file=adjustl(trim(filein))//trim(adjustl('.RTspop.fin')),status='old')
do i=1,nheader
  read(2,*) skipline
enddo

!continue reading RTspop.fin
do p=1,itot
  do j=1,coo
      read(2,'(100ES15.7)',advance='no')  pdr(p)%coolant(j)%pop(1:coolant(j)%cnlev)
  enddo
  read(2,*)
enddo
close(2)

return
end subroutine

