program init3DPDR
use m_Mesh
use m_Ray_box
implicit none
real*8,parameter:: Pi = 3.141592653589793238462643383279502884
integer :: i,j,k,p,nunit
real*8 :: x,y,z,xc,yc,zc,r,rr,n
real(8),allocatable,dimension(:,:,:) :: density_3D
real(8),allocatable,dimension(:) :: density_1D

! locals
real(RK) :: xs,ys,zs
integer :: levels
type(box) :: box1
type(HEALPix_ray) :: ray
real(RK) :: thfpix, phfpix, contribution, contribution_1D
real(RK) :: corner_min(3), corner_max(3)

xlx = 3.D0
yly = 3.D0 
zlz = 3.D0
nxc = 32
nyc = 32
nzc = 32
dx = xlx/real(nxc,kind=RK)
dy = yly/real(nyc,kind=RK)
dz = zlz/real(nzc,kind=RK)

xc = 1.5D0
yc = 1.5D0
zc = 1.5D0
r = 1.D0

allocate(density_3D(nxc,nyc,nzc))
allocate(density_1D(nxc*nyc*nzc))

open(nunit,file='3Dsphere.dat')
do k = 1,nzc
  do j = 1,nyc
    do i = 1,nxc
        x = (DBLE(i)-0.5D0)*dx
        y = (DBLE(j)-0.5D0)*dy
        z = (DBLE(k)-0.5D0)*dz
        rr = sqrt((x-xc)**2.D0+(y-yc)**2.D0+(z-zc)**2.D0)
        if(rr.le.r) then
          n=100.D0
        else
          n=0.D0
        endif
        write(nUnit,*) x,y,z,n

        p = I + (J-1)*nxc + (K-1)*(nxc*nyc) 
        density_3D(i,j,k) = n
        density_1D(p) = n 
    enddo
  enddo
enddo
close(nUnit)

levels = nint(log(DBLE(nxc)) / log(2.D0)) + 1
corner_min = [0.D0, 0.D0, 0.D0]
corner_max = [xlx, yly, zlz]
box1%min = corner_min
box1%max = corner_max
xs = 1.5D0 - 0.5D0*dx
ys = 1.5D0 - 0.5D0*dy
zs = 1.5D0 - 0.5D0*dz
thfpix = 1.5707963267948966D0
phfpix = 4.7123889803846897D0

ray%origin = [xs, ys, zs]
ray%angle = [thfpix, phfpix]
print*, thfpix*180.D0/pi, phfpix*180.D0/pi

contribution = 0.D0
contribution_1D = 0.D0
call octree(ray, box1, levels-1, contribution, contribution_1D,density_3D,density_1D)
print*,'3D',contribution
print*,'1D',contribution_1D
print*,'size',dx

end program
