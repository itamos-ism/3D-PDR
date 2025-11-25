subroutine writeparams
  use definitions
  use healpix_types
  use maincode_module, only : input, level, Tguess, iterstep,&
                              & theta_crit, ITERTOT, &
                              & Gext, AV_fac, UV_fac, nspec, nreac, maxpoints, &
                              & Tlow0, Thigh0, Tmin, Tmax, Fcrit, Tdiff, dust_temperature,&
                              & chemiterations, zeta,  v_turb!, fieldchoice
  use global_module, only : metallicity, omega, grain_radius


write(6,*) 'Input file:               ',input
write(6,*) 'HEALPix level:            ',level
#ifdef ONEDIMENSIONAL
if (level.gt.0) STOP "HEALPix level must be set to 0 in ONEDIMENSIONAL mode"
#endif
#ifdef PSEUDO_2D
if (level.gt.0) STOP "HEALPix level must be set to 0 in ONEDIMENSIONAL mode"
#endif
write(6,*) 'Theta critical:           ',theta_crit
write(6,*) 'Angle between rays:       ',sqrt(pi/3.0D0/4.0D0**(real(level)))
write(6,*) 'Maxpoints                 ',maxpoints
if (theta_crit.ge.pi/2.0D0) stop 'theta_crit must be less than pi/2'
write(6,*) 'Guess Temperature (K):    ',Tguess
write(6,*) 'Dust  Temperature (K):    ',dust_temperature
write(6,*) 'Turbulent velocity (cm/s):',v_turb
#ifdef THERMALBALANCE
write(6,*) 'Tlow:                     ',Tlow0
write(6,*) 'Thigh:                    ',Thigh0
#endif
write(6,*) 'Tmin:                     ',Tmin
write(6,*) 'Tmax:                     ',Tmax
write(6,*) 'Fcrit:                    ',Fcrit
write(6,*) 'Tdiff:                    ',Tdiff
!write(6,*) 'Form of field:            ',fieldchoice
write(6,*) 'Gext:                     ',Gext(1)
write(6,*) 'AV factor:                ',AV_fac
write(6,*) 'UV factor:                ',UV_fac
#ifdef REDUCED
write(6,*) 'Chemical network:           REDUCED'
#elif MEDIUM
write(6,*) 'Chemical network:           MEDIUM'
#elif FULL
write(6,*) 'Chemical network:           FULL'
#elif MYNETWORK
write(6,*) 'Chemical network:           MYNETWORK'
#endif
write(6,*) 'Number of species:        ',nspec
write(6,*) 'Number of reactions:      ',nreac
write(6,*) 'Total iterations:         ',itertot
write(6,*) 'Output interval / iter.:  ',iterstep
write(6,*) 'Chemiterations:           ',chemiterations
write(6,*) 'Zeta:                     ',zeta*1.3d-17
write(6,*) 'Metallicity               ',metallicity
write(6,*) 'Omega                     ',omega
write(6,*) 'Grain radius              ',grain_radius
close(12)
write(6,*) '============================================='
write(6,*) ''
write(6,*) '------FLAGS-----'
#ifdef THERMALBALANCE
write(6,*) 'THERMALBALANCE'
#endif
#ifdef GRAINRECOMB
write(6,*) 'GRAINRECOMB'
#endif
#ifdef SUPRATHERMAL
write(6,*) 'SUPRATHERMAL'
#endif
#ifdef OPENMP
write(6,*) 'OPENMP'
#endif
#ifdef ONEDIMENSIONAL
write(6,*) 'ONEDIMENSIONAL'
#else
write(6,*) 'FULL_3D'
#endif
#ifdef DUST
write(6,*) 'DUST'
#endif
#ifdef FORCECONVERGENCE
write(6,*) 'FORCECONVERGENCE'
#endif
#ifdef H2FORM
write(6,*) 'H2FORM'
#endif
#ifdef GUESS_TEMP
write(6,*) 'GUESS_TEMP'
#endif
write(6,*) ''
write(6,*) 'Reading initial conditions file'

return
end subroutine
