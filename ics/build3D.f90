program one
implicit none
integer::i,j,k,ctot
real::x,y,z,r,n,sbox,csize,sradius,rho

write(6,*) 'Grid resolution per axis (e.g. 32 for 32^3)'
read(5,*) ctot
write(6,*) 'Box size in units of pc (e.g. 10)'
read(5,*) sbox
write(6,*) 'Radius of sphere in pc (e.g. 5)'
read(5,*) sradius
write(6,*) 'Density of sphere in cm-3 (e.g. 1000)'
read(5,*) rho

!ctot = 32               !grid resolution
!sbox = 10.               !grid size
!sradius = 5.            !radius of sphere
!rho = 1000.             !density of sphere
csize = sbox/real(ctot) !cell size

open(unit=12,file='3Dsphere.dat',status='replace')

write(12,*) ctot, ctot, ctot
write(12,*) sbox, sbox, sbox

do i=1,ctot
  x = sbox*real(i)/real(ctot)-csize/2.-sbox/2.
  do j=1,ctot
    y = sbox*real(j)/real(ctot)-csize/2.-sbox/2.
    do k=1,ctot
      z = sbox*real(k)/real(ctot)-csize/2.-sbox/2.

      r = sqrt(x**2 + y**2 + z**2)

      n=0.0
      if (r.le.sradius) n = rho

      write(12,*) x+sbox/2.,y+sbox/2.,z+sbox/2.,n

    enddo
  enddo
enddo

end program

