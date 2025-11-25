program convert_fin2h5
use m_paramters
use m_visual
use m_writeoutputs
implicit none

integer :: pdr_ptot,unit
integer :: i,j,k,II,JJ,KK,p,pp,ilevel,ico
character(len=7) :: first_col
logical :: file_exists

pdr_ptot = nxc*nyc*nzc

allocate(pdr(nxc,nyc,nzc))
do k=1,nzc
do j=1,nyc
do i=1,nxc
  allocate(pdr(i,j,k)%coolant(1:coo))
  allocate(pdr(i,j,k)%AV(0:nrays-1))
  allocate(pdr(i,j,k)%abundance(1:nspec))
  do ii=1,coo
    allocate(pdr(i,j,k)%cooling(ii))
    allocate(pdr(i,j,k)%coolant(ii)%pop(cnlev(ii)))
  enddo
  do ii=1,12
    allocate(pdr(i,j,k)%heating(ii))
  enddo
enddo
enddo
enddo

! read fin
fin_file = "../"//trim(adjustl(outname))//".pdr.fin"
inquire(file=trim(fin_file), exist=file_exists)
if (file_exists) then
  open(newunit=unit,file=fin_file,status='old')
  do i=1,nxc
  do j=1,nyc
  do k=1,nzc
      read(unit,'(A7,5ES15.7,I5,300ES15.7)') first_col,pdr(i,j,k)%x, pdr(i,j,k)%y, pdr(i,j,k)%z, pdr(i,j,k)%Tgas,pdr(i,j,k)%Tdust,&
      &pdr(i,j,k)%etype,pdr(i,j,k)%rho,pdr(i,j,k)%UVfield,pdr(i,j,k)%abundance,pdr(i,j,k)%AV
  enddo
  enddo
  enddo
  close(unit)
endif

fin_file = "../"//trim(adjustl(outname))//".cr.fin"
inquire(file=trim(fin_file), exist=file_exists)
if (file_exists) then
  open(newunit=unit,file=fin_file,status='old')
  do i=1,nxc
  do j=1,nyc
  do k=1,nzc
      read(unit,'(ES15.7)') pdr(i,j,k)%zetalocal
  enddo
  enddo
  enddo
  close(unit)
endif

fin_file = "../"//trim(adjustl(outname))//trim(adjustl(".cool"))//".fin"
inquire(file=trim(fin_file), exist=file_exists)
if (file_exists) then
  open(newunit=unit,file=fin_file,status='old')
  do i=1,nxc
  do j=1,nyc
  do k=1,nzc
      read(unit,'(A7,200ES15.7)') first_col, pdr(i,j,k)%x, pdr(i,j,k)%y, pdr(i,j,k)%z, pdr(i,j,k)%cooling(:),pdr(i,j,k)%totalcooling, pdr(i,j,k)%AV(:)
  enddo
  enddo
  enddo
  close(unit)
endif

fin_file = "../"//trim(adjustl(outname))//trim(adjustl(".heat"))//".fin"
inquire(file=trim(fin_file), exist=file_exists)
  if (file_exists) then
  open(newunit=unit,file=fin_file,status='old')
  do i=1,nxc
  do j=1,nyc
  do k=1,nzc
      read(unit,'(A7,200ES15.7)') first_col, pdr(i,j,k)%x, pdr(i,j,k)%y, pdr(i,j,k)%z, pdr(i,j,k)%heating, pdr(i,j,k)%AV
  enddo
  enddo
  enddo
  close(unit)
endif

fin_file = "../"//trim(adjustl(outname))//trim(adjustl(".RTspop"))//".fin"
inquire(file=trim(fin_file), exist=file_exists)
if (file_exists) then
  open(newunit=unit,file=fin_file,status='old')
  ico=0
  do
    read(unit,'(A)') cfile
    cfile=trim(adjustl(cfile))
    if (cfile=='ENDCOOLFILES') exit
    ico=ico+1
    coolfile(ico)=cfile
  enddo
  ! print*,coolfile
  
  do i=1,nxc
  do j=1,nyc
  do k=1,nzc
      do ico = 1, coo
          read(unit, '(100ES15.7)', advance='no') pdr(i,j,k)%coolant(ico)%pop(1:cnlev(ico))
      end do
      read(unit, *)
  enddo
  enddo
  enddo
  close(unit)
endif

fin_file = "../"//trim(adjustl(outname))//trim(adjustl(".velocity"))//".dat"
inquire(file=trim(fin_file), exist=file_exists)
if (file_exists) then
  open(newunit=unit,file=fin_file,status='old')
  do i=1,nxc
  do j=1,nyc
  do k=1,nzc
      read(unit,*) pdr(i,j,k)%vx, pdr(i,j,k)%vy, pdr(i,j,k)%vz
  enddo
  enddo
  enddo
  close(unit)
endif
print*, 'Finish reading fin'

! write h5
call writeoutputs
print*, 'Finish writing h5'

! write xmf
dx = xlx/real(nxc,kind=RK)
dy = yly/real(nyc,kind=RK)
dz = zlz/real(nzc,kind=RK)
call initvisual

! print*,pdr(nxc,nyc,nzc-1)%coolant(4)%pop(cnlev(4)-1)
print*, 'Done'

end program convert_fin2h5
