module m_paramters
  implicit none

  integer, parameter :: RK = KIND(0.D0)
  integer, parameter :: coo = 4
  integer, parameter :: nspec = 33
  integer, parameter :: nrays = 12
  integer, parameter, dimension(coo) :: cnlev = [41, 5, 5, 5]
  integer, parameter :: nxc = 32, nyc = 32, nzc = 32
  real(RK), parameter :: xlx = 3.94156170, yly = 3.93756294, zlz = 3.93856335
  character(len=50), parameter :: outname = 'test'
  
  real(RK) :: dx,dy,dz
  character(len=20) :: cfile, coolfile(1:10)

  type pdr_excit
    real(RK), pointer :: pop(:)
  end type pdr_excit

  type pdr_node
      integer :: etype
      real(RK) :: x,y,z,vx,vy,vz,Tgas,Tdust,rho,UVfield,totalcooling
      real(RK), pointer :: AV(:)
      real(RK), pointer :: abundance(:)
      real(RK), pointer :: cooling(:)
      real(RK), pointer :: heating(:)
      type(pdr_excit),allocatable :: coolant(:)
      real(RK) :: zetalocal
  end type pdr_node
  type(pdr_node), allocatable :: pdr(:,:,:)

end module m_paramters
