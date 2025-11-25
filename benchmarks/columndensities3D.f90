program one
implicit none
integer::i,j,k,id,etype,ll,itot
double precision::pc2cm,mp,N(1:33),Tgas,Tdust,uv
double precision,allocatable::nH(:,:,:),x(:),y(:),z(:)
double precision,allocatable::abun(:,:,:,:)
character(len=50)::prefix,input,output
integer :: spec, sCO, sC, sCp, sH, sH2
character(len=1)::net

write(6,*) 'Give grid resolution (e.g. 64)'
read(5,*) itot
write(6,*) 'Which network did you use? (R: reduced, M: medium, F: full)'
read(5,*) net
if (net.eq.'R') then
      write(6,*) 'Reduced network selected'
      spec = 33
      sCO = 28
      sC  = 25
      sCp = 11
      sH  = 32
      sH2 = 31
else if (net .eq. 'M') then 
      write(6,*) 'Medium network selected'
      spec = 77
      sCO = 28
      sC  = 8
      sCp = 14
      sH  = 9
      sH2 = 7
else if (net .eq. 'F') then
      write(6,*) 'Full network selected'
      spec = 215
      sCO = 211
      sC  = 210
      sCp = 172
      sH  = 214
      sH2 = 213
else 
      STOP 'Invalid network; try again'
end if


allocate(nH(1:itot,1:itot,1:itot))
allocate(abun(1:spec,1:itot,1:itot,1:itot))
allocate(x(1:itot),y(1:itot),z(1:itot))
write(6,*) 'give prefix'
read(5,*) prefix
input = trim(adjustl(prefix))//".pdr.fin"
output = 'cds_'//trim(adjustl(prefix))//".dat"
open(unit=1,file=input,status='old')
write(6,*) 'Reading input'
do i=1,itot
  write(6,*) i
  do j=1,itot
    do k=1,itot
      read(1,*) id,x(i),y(j),z(k),Tgas,Tdust,etype,nH(i,j,k),uv,abun(1:spec,i,j,k)
    enddo
  enddo
enddo
close(1)

pc2cm = 3.0856d18

ll=0

open(unit=2,file=output,status='replace')
write(6,*) 'Writing output'
do i=1+ll,itot-ll
  write(6,*) i
  do j=1+ll,itot-ll
    N = 0
    do k=1+ll,itot-1-ll
       N = N + 0.5*(nH(i,j,k)*abun(:,i,j,k)+nH(i,j,k+1)*abun(:,i,j,k+1))*abs(z(k+1)-z(k))*pc2cm
    enddo
    write(2,*) x(i),y(j),N(sH2),N(sH),N(sCp),N(sC),N(sCO)
  enddo
  write(2,*) ''
enddo
close(2)


end program
