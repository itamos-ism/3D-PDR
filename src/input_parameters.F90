
SUBROUTINE readparams
 
!T.Bisbas, T.Bell
  use definitions
  use healpix_types
  use maincode_module, only : input, level, Tguess, v_turb_inp, &
                              & theta_crit, ITERTOT, output, &
                              & Gext, AV_fac, UV_fac, nspec, nreac, maxpoints,&
                              & Tlow0, Thigh0, Tmin, Tmax, Fcrit, Tdiff, dust_temperature,&
                              & chemiterations, zeta, end_time, Av_crit, v_alfv, &
                              & indir, outdir, coolfile, coo, redshift, Tcmb, crfieldchoice, &
                              & paramFile, crattennorm, crattenslope, crattenn0, UVdirchoice, &
                              & user_UVAngle, UVdir
  use global_module, only : metallicity, omega, grain_radius
  use chemistry_module
  use m_Mesh


  real(kind=dp):: dummy
  character(len=50)::coolin,suffix

open(unit=12,file=paramFile,status='old')
read(12,*); read(12,*); read(12,*)
read(12,*) indir
read(12,*) input
read(12,*) outdir
read(12,*) output
input = trim(adjustl(indir))//'/'//trim(adjustl(input))
output = trim(adjustl(outdir))//'/'//trim(adjustl(output))
#ifdef RAYTHEIA
open(unit=2,file=input,status='old')
read(2,*) nxc, nyc, nzc
read(2,*) xlx, yly, zlz
close(2)
#endif
read(12,*); read(12,*); read(12,*)
Gext=0
read(12,*) Gext(1)
#ifdef CRATTENUATION
#if CRATTENUATION == 1
read(12,*) crfieldchoice !L or H
select case (crfieldchoice)
  case ('L', 'l', 'H', 'h')
      continue
  case default
      STOP 'Invalid crfieldchoice'
end select
#elif CRATTENUATION == 2
read(12,*) crattennorm
read(12,*) crattenn0
read(12,*) crattenslope
#else
  write(6,*) "Invalid CRATTENUATION Flag. Must be 1 or 2"
  stop
#endif
#else
read(12,*) zeta
zeta = zeta/1.3d-17
#endif
read(12,*) metallicity
read(12,*) v_turb_inp
read(12,*) end_time
read(12,*) grain_radius
omega=0.42
read(12,*) Av_fac
AV_fac = AV_fac*metallicity
read(12,*) UV_fac
read(12,*) redshift
Tcmb = 2.725*(1.+redshift)
read(12,*) Av_crit
read(12,*) v_alfv
v_alfv = (v_alfv**2) * 40.382518!K --> protonmass * (1 km/s)^2 / 3 / k_boltzmann
read(12,*); read(12,*); read(12,*)
read(12,*) level
read(12,*) theta_crit
read(12,*); read(12,*); read(12,*)
read(12,*) relative_abundance_tolerance
read(12,*) absolute_abundance_tolerance
read(12,*); read(12,*); read(12,*)
read(12,*) chemiterations
read(12,*) ITERTOT
read(12,*); read(12,*); read(12,*)
read(12,*) Tguess
read(12,*) Tmin
if (Tmin<Tcmb) Tmin=Tcmb
read(12,*) Tmax
read(12,*) dust_temperature
read(12,*) Fcrit
read(12,*) Tdiff
Tlow0=10.0
Thigh0=8000.0
read(12,*); read(12,*); read(12,*)
coo=0
do
  read(12,*,end=99) coolin
  coo=coo+1
  coolfile(coo) = coolin
enddo
99 continue
coolfile = 'chemfiles/'//coolfile

#ifdef REDUCED
suffix = 'reduced.d'
#elif MEDIUM
suffix = 'medium.d'
#elif FULL
suffix = 'full.d'
#elif MYNETWORK
suffix = 'mynetwork.d'
#endif
close(1)
open(unit=1,file='chemfiles/species_'//trim(adjustl(suffix)),status='old')
nspec=0
do 
 read(1,*,end=100) dummy
 nspec=nspec+1
enddo
100 continue
close(1)
open(unit=1,file='chemfiles/rates_'//trim(adjustl(suffix)),status='old')
nreac=0
do 
 read(1,*,end=101) dummy
 nreac=nreac+1
enddo
101 continue
close(1)

#ifdef RAYTHEIA
call Mesh
#endif

close(12)
return
END SUBROUTINE readparams
