module RTtool_module

double precision, parameter, public ::  KB = 1.38065040D-16 !Boltzmann constant cgs
double precision, parameter, public ::  C  = 2.99792458D+10 !speed of light cgs
double precision, parameter, public ::  MP = 1.67262164D-24 !proton mass cgs
double precision, parameter, public ::  HP = 6.62606896D-27 !Planck's constant cgs
double precision, parameter, public ::  HB = 1.05457163D-27 !Planck's constant / 2pi
double precision, parameter, public ::  HK = 4.79923734D-11 !Planck's constant / Boltzmann constant
double precision, parameter, public ::  NA = 6.02214179D+23 !Avogadro's number
double precision, parameter, public ::  AU = 1.66053878D-24 !atomic mass unit
double precision, parameter, public ::  MH = 1.67372346D-24 !hydrogen mass cgs
double precision, parameter, public ::  ME = 9.10938215D-28 !electron mass cgs
double precision, parameter, public ::  EC = 4.80320427D-10 !elementary charge in esu
double precision, parameter, public ::  PC = 3.08568025D+18 !pc in cm
double precision, parameter, public ::  EV = 1.60217646D-12 !electron volt in erg
double precision, parameter, public ::  PI = 3.141592653589

real :: Tcmb, z

integer :: il,jl
integer :: iH2
integer :: coolev
character(len=20) :: suffix
character(len=20) :: skipline
integer :: nheader
#ifdef LEVMAX
integer :: lmax_levels = 0
#endif
integer::p,i,coo
character(len=30) :: coolfile(1:30)
character(len=40) :: network
character(len=50) :: pdrfile,popfile
type coolant_node
  double precision, pointer :: COEFF(:)
  double precision, pointer :: ENERGIES(:), WEIGHTS(:)
  double precision, pointer :: A_COEFFS(:,:), B_COEFFS(:,:), C_COEFFS(:,:)
  double precision, pointer :: FREQUENCIES(:,:), TEMPERATURES(:,:)
  double precision, pointer :: HP_COL(:,:,:)
  double precision, pointer :: H_COL(:,:,:)
  double precision, pointer :: EL_COL(:,:,:)
  double precision, pointer :: HE_COL(:,:,:)
  double precision, pointer :: H2_COL(:,:,:)
  double precision, pointer :: PH2_COL(:,:,:)
  double precision, pointer :: OH2_COL(:,:,:)
  double precision :: molweight
  double precision :: N
  integer :: cnlev, cntemp, cspec
  character(len=10) :: cname
  double precision, pointer :: Tr(:,:)
  double precision, pointer :: units(:,:) !unit conversion to K
end type coolant_node
type(coolant_node), allocatable::coolant(:)
integer::cur_nlev,cur_ntemp

character(len=40)::chemsuf
character(len=50)::prefix,directory,filein
character(len=10),allocatable::species(:)

integer::id,etype,itot,nspec
real::uv,dummy

double precision :: Ntot
double precision, allocatable :: Nabn(:)

type pdr_excit
   double precision, pointer :: pop(:)            
   double precision, pointer :: Tex(:,:), Bnu(:,:), tau(:,:)
   double precision, pointer :: Tr(:,:)
   double precision :: N
end type pdr_excit

type pdr_node
   type(pdr_excit), allocatable :: coolant(:)
   double precision, pointer :: abun(:)      
   double precision, pointer :: Nabn(:)
   double precision :: Av
   double precision :: Ntot
   double precision :: rho    
   double precision :: x      
   double precision :: Tgas,Tdust
end type pdr_node
type (pdr_node), allocatable :: pdr(:)      !main 3DPDR array for each grid point p

double precision :: step,sigma,phi,frac
double precision :: dtau, Tr_incr1, Tr_incr2
double precision :: vturb
double precision :: avg_nH2, avg_Tgas, Ntgas, Nnh2
double precision :: Bckgr


end module RTtool_module
