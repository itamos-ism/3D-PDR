subroutine readdensity
use definitions
use healpix_types
use healpix_module
use maincode_module
use global_module
#ifdef OPENMP
use omp_lib
#endif
use chemistry_module
use maincode_local
use m_Mesh
use m_Ray_box

integer :: temp


open(unit=2,file=input,status='old')

#ifndef ONEDIMENSIONAL
read(2,'(/)')
#endif
#ifdef RAYTHEIA
pdr_ptot = nxc*nyc*nzc
#else
pdr_ptot = 0
do 
  read(2,*,end=100) xpos
  pdr_ptot = pdr_ptot + 1
enddo
100 continue
rewind(2)
#endif
allocate(pdr(1:pdr_ptot))
do p=1,pdr_ptot
    read(2,*) xpos,ypos,zpos,denst
    pdr(p)%etype = 1
    pdr(p)%x=xpos
    pdr(p)%y=ypos
    pdr(p)%z=zpos
    pdr(p)%rho=denst
enddo
write(6,*) 'PDR elements        = ',pdr_ptot
write(6,*) 'Maximum PDR density = ',maxval(pdr(:)%rho)
write(6,*) 'Minimum PDR density = ',minval(pdr(:)%rho)
#ifdef RAYTHEIA_MO
write(6,*) 'Grid resolution     = ',nxc
write(6,*) 'Domain size         = ',real(xlx)
nside=2**level
nrays=12*nside**2
ns_max=8192
levels = nint(log(DBLE(nxc)) / log(2.D0)) + 1
corner_min = [0.D0, 0.D0, 0.D0]
corner_max = [xlx, yly, zlz]

allocate(maxpoints_ray(0:nrays-1))
maxpoints = 0
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p, j, box1, thfpix, phfpix, ray, maxpoints_ray) REDUCTION(MAX: maxpoints)
#endif
do p=1,pdr_ptot
  box1%min = corner_min
  box1%max = corner_max
  maxpoints_ray = 0
  do j=0,nrays-1
      call pix2ang_nest(nside, j, thfpix, phfpix)
      ray%origin = [pdr(p)%x, pdr(p)%y, pdr(p)%z]
      ray%angle = [thfpix, phfpix]

      call raytheia_maxpoints(ray, box1, levels-1, j, maxpoints_ray)
  enddo
  temp = maxval(maxpoints_ray)
  if(temp.gt.maxpoints) maxpoints = temp 
enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
print*,'Maxpoints           = ',maxpoints

#else
#ifdef RAYTHEIA
write(6,*) 'Grid resolution     = ',nxc
write(6,*) 'Domain size         = ',real(xlx)
levels = nint(log(DBLE(nxc)) / log(2.D0)) + 1
corner_min = [0.D0, 0.D0, 0.D0]
corner_max = [xlx, yly, zlz]
maxpoints = 3*nxc - 2 
#else
maxpoints = 1000
allocate(rra(0:pdr_ptot))
allocate(rrb(1:pdr_ptot))
allocate(pdrpoint(1:3,1:pdr_ptot)) 
close(2)
do p=1,pdr_ptot
  rra(p) = sqrt(pdr(p)%x**2+pdr(p)%y**2+pdr(p)%z**2)
  rrb(p) = p
  pdrpoint(1,p) = pdr(p)%x
  pdrpoint(2,p) = pdr(p)%y
  pdrpoint(3,p) = pdr(p)%z
enddo
#endif
#endif

return
end subroutine
