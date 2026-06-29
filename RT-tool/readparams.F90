subroutine readparams
use RTtool_module
character(len=200)::cfile
character(len=10)::cname
integer::cindex

close(1);open(unit=1,file='paramsRT.dat',status='old')

read(1,*) directory
read(1,*) prefix
filein = adjustl(trim(directory))//'/'//adjustl(trim(prefix))

! The chemical network is read from src/config.mk (the current build) rather than
! from the .RTspop.fin header: the header line is missing for custom networks
! (e.g. ones built with the Network Builder), so reading it there is unreliable.
! Only the species-file suffix is derived per network; nspec and the H2 index are
! read straight from the species file so RT-tool always matches the network.
call rt_read_network(network)
call rt_chemsuf(network, chemsuf)

close(3); open(unit=3,file='chemfiles/species_'//trim(adjustl(chemsuf)),status='old')
nspec = 0
iH2 = 0
do
  read(3,*,end=50) cindex, cname
  nspec = nspec + 1
  if (trim(adjustl(cname)).eq.'H2') iH2 = nspec
enddo
50 continue
close(3)
if (iH2.eq.0) stop 'ERROR: H2 species not found in species file'
write(6,*) 'Number of species = ',nspec,' (H2 index = ',iH2,')'

allocate(species(1:nspec))
write(6,*) 'Chemical network = ',trim(adjustl(network))

read(1,*) vturb
vturb = vturb*1d5 !km/s to cm/s
read(1,*) z
Tcmb = 2.725*(1+z)

! Read the coolant files from the .RTspop.fin header. Any leading line that is
! not a coolant file (e.g. an optional network header line) is skipped, so this
! works whether or not the header carries a network line. nheader counts ALL
! header lines (up to and including ENDCOOLFILES) so the data are read correctly.
open(unit=2,file=adjustl(trim(filein))//trim(adjustl('.RTspop.fin')),status='old')
nheader = 0
coo = 0
do
  read(2,'(A)') cfile
  cfile = trim(adjustl(cfile))
  nheader = nheader + 1
  if (cfile=='ENDCOOLFILES') exit
  if (index(cfile,'.dat') > 0) then
    coo = coo + 1
    coolfile(coo) = cfile
    write(6,*) trim(coolfile(coo))
  endif
#ifdef LEVMAX
  if (index(cfile,'LMAX=') == 1) then
    read(cfile(6:),*) lmax_levels
    write(6,'(a,i4)') ' [LEVMAX] RTspop: level cap from header: ', lmax_levels
  endif
#endif
enddo
close(2)
write(6,*) 'Coolants found = ',coo

return
end subroutine


! ---------------------------------------------------------------------------
! Read NETWORK (and XRAYS) from src/config.mk and return the network name.
! ---------------------------------------------------------------------------
subroutine rt_read_network(net)
character(len=*),intent(out)::net
character(len=300)::line
character(len=40)::xrays
integer::ios,eq
net = ''
xrays = '0'
open(unit=9,file='src/config.mk',status='old',iostat=ios)
if (ios/=0) open(unit=9,file='../src/config.mk',status='old',iostat=ios)
if (ios/=0) return
do
  read(9,'(A)',iostat=ios) line
  if (ios/=0) exit
  line = adjustl(line)
  if (len_trim(line)==0) cycle
  if (line(1:1)=='#') cycle
  eq = index(line,'=')
  if (eq<=0) cycle
  if (index(line(1:eq),'NETWORK') > 0) net   = trim(adjustl(line(eq+1:)))
  if (index(line(1:eq),'XRAYS')   > 0) xrays = trim(adjustl(line(eq+1:)))
enddo
close(9)
! X-ray reduced network uses a different species file
if (trim(net)=='REDUCED' .and. trim(xrays)=='1') net = 'Xreduced'
end subroutine


! ---------------------------------------------------------------------------
! Map a network name to its species-file suffix (species_<suffix>).
! Built-in networks use their lowercase file names; any other (custom) network
! uses its own name verbatim, matching the Network Builder's species_<NAME>.d.
! ---------------------------------------------------------------------------
subroutine rt_chemsuf(net, suf)
character(len=*),intent(in)::net
character(len=*),intent(out)::suf
select case (trim(adjustl(net)))
case ('REDUCED');   suf = 'reduced.d'
case ('Xreduced');  suf = 'Xreduced.d'
case ('MEDIUM');    suf = 'medium.d'
case ('MEDIUMS');   suf = 'mediumS.d'
case ('SCUTUM');    suf = 'scutum.d'
case ('FULL');      suf = 'full.d'
case ('MYNETWORK'); suf = 'mynetwork.d'
case default;       suf = trim(adjustl(net))//'.d'
end select
end subroutine
