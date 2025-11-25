program logarithmic

real::avmax,avmin
real::lambda,n
real::x,av,y,z,size,yzdim
integer::i,j,tot,resolution
real::av0,cm,length
!************** V
!real::av_1,av_2,displacement
!integer::NN
!************** ^

cm=3.0856e18
av0=6.289e-22
write(6,*) 'give number density [cm^-3] (e.g. 1000)'
read(5,*) n
write(6,*) 'give avmax (e.g. 10)'
read(5,*) avmax
length=avmax/cm/av0/n
write(6,*) 'length [pc] = ',length
write(6,*) 'give log10(avmin) (e.g. -3)'
read(5,*) avmin
10 continue
write(6,*) 'give x-resolution (e.g. 30)'
read(5,*) resolution

tot=0
open(unit=1,file='outgrid.dat',status='replace')
write(1,'(4E15.7)') 0.,0.,0.,n

lambda=(log10(avmax)-avmin)*real(resolution)
do i=0,int(lambda)
  av=10.0**(avmin+real(i)/real(resolution))
  x=av/cm/av0/n
  write(1,'(4E15.7)') x,0.,0.,n
  tot=tot+1
enddo

write(6,*) 'file [outgrid.dat] written!'
end program
