module m_Mesh
  implicit none
  private
  integer,parameter,public::RK=KIND(0.D0)

  ! mesh option
  real(RK),public::xlx,yly,zlz ! domain length
  integer,public::nxp,nyp,nzp ! grid points number
  integer,public::nxc,nyc,nzc ! grid center number
  real(RK),public::dx,dy,dz

  public::Mesh
contains

  subroutine Mesh
    implicit none

    ! locals
    integer:: nUnit,ierror

    ! grid parameters
    nxp=nxc+1
    nyp=nyc+1
    nzp=nzc+1
    dx = xlx/real(nxc,kind=RK)
    dy = yly/real(nyc,kind=RK)
    dz = zlz/real(nzc,kind=RK)
  end subroutine Mesh

end module m_Mesh
