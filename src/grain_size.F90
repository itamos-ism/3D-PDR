#ifdef MRNDUST
!=======================================================================
!
!  grain_distribution_module - MRN (Mathis, Rumpl & Nordsieck 1977) grain
!  size distribution, used by the dust-size-dependent rates of the code
!  when MRNDUST = 1 (GRAINRECOMB = 2, H2 formation, gas-grain heating).
!
!  dn/da = A_MRN * a**mrn_q      for mrn_amin <= a <= mrn_amax
!
!  A_MRN is fixed, at each grid point, by requiring the total grain MASS
!  density to reproduce the model's dust-to-gas mass ratio
!  (1.0D-2*metallicity, the "Dust-to-gas normalized to 1e-2" params.dat
!  entry) at the local gas density:
!
!    rho_dust = INT[ (4/3) pi a^3 rho_grain * dn/da da ]
!             = (1.0D-2*metallicity) * density * MH
!
!  All routines are pure functions without module-level mutable state
!  (the log-spaced size grid of MRN_DS_G is a local array rebuilt at each
!  call), so they can be called safely from inside OpenMP parallel regions.
!
!-----------------------------------------------------------------------
MODULE grain_distribution_module

  USE healpix_types, ONLY : dp, i4b, PI, KB, MH, EC

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: MRN_SIGMA_PER_NH, MRN_DS_G, &
          & mrn_amin, mrn_amax, mrn_q, mrn_nbin, mrn_rhogr

! -------- MRN distribution parameters --------
! Standard MRN (1977) values: a_min and a_max span 5 nm - 0.25 micron
! (no PAHs). mrn_rhogr is the grain material density (3 g cm^-3, the same
! value used by the single-grain GRAINRECOMB = 2 treatment).
  REAL(kind=dp), PARAMETER :: mrn_amin  = 5.0D-7   ! cm, min grain radius (5 nm)
  REAL(kind=dp), PARAMETER :: mrn_amax  = 2.5D-5   ! cm, max grain radius (0.25 um)
  REAL(kind=dp), PARAMETER :: mrn_q     = -3.5D0   ! MRN power-law index
  INTEGER(kind=i4b), PARAMETER :: mrn_nbin = 40    ! log-spaced size bins (DS87 integral)
  REAL(kind=dp), PARAMETER :: mrn_rhogr = 3.0D0    ! g cm^-3, grain material density

CONTAINS

!-----------------------------------------------------------------------
!  Total grain geometric cross-section area per H nucleus (cm^2),
!  i.e. sum over grains of pi a^2 per H nucleus, the MRN-integrated
!  replacement for the hard-wired per-H-nucleus grain cross sections.
!  Density independent by construction (A_MRN ~ density cancels against
!  the explicit density in n_gr*sigma_gr/density).
!-----------------------------------------------------------------------
  FUNCTION MRN_SIGMA_PER_NH(metallicity) RESULT(sigma_per_nh)
    REAL(kind=dp), INTENT(IN) :: metallicity
    REAL(kind=dp) :: sigma_per_nh
    REAL(kind=dp) :: normint_mass, normint_area

    normint_mass = (mrn_amax**(mrn_q+4.0D0) - mrn_amin**(mrn_q+4.0D0)) / (mrn_q+4.0D0)
    normint_area = (mrn_amax**(mrn_q+3.0D0) - mrn_amin**(mrn_q+3.0D0)) / (mrn_q+3.0D0)

    sigma_per_nh = 0.75D0*(1.0D-2*metallicity)*MH*normint_area / (mrn_rhogr*normint_mass)

  END FUNCTION MRN_SIGMA_PER_NH

!-----------------------------------------------------------------------
!  Draine & Sutin (1987) neutral-grain (polarization-limit) geometric
!  rate factor, integrated over the MRN size grid instead of evaluated at
!  one representative grain radius:
!
!    G(T,Z) = SUM_bins  n_gr(a) * pi*a^2 * J(tau(a,Z)) * da   [cm^-1]
!    n_gr(a) da = A_MRN * a**mrn_q * da
!    tau(a,Z) = a*KB*T / (Z*EC)^2
!    J(tau) = 1 + sqrt(pi/(2*tau))
!
!  so that k_gr(X^Z+) = G(T,Z) * v_X (thermal speed of the ion), as in the
!  single-grain treatment. The bin sum (midpoint rule on mrn_nbin
!  log-spaced bins) is needed because J(tau(a)) is nonlinear in a.
!-----------------------------------------------------------------------
  FUNCTION MRN_DS_G(metallicity, density, temperature, zion) RESULT(G)
    REAL(kind=dp), INTENT(IN) :: metallicity, density, temperature
    INTEGER(kind=i4b), INTENT(IN) :: zion
    REAL(kind=dp) :: G

    REAL(kind=dp) :: normint_mass, A_mrn
    REAL(kind=dp) :: bin_a, bin_da, tau, Jt
    INTEGER(kind=i4b) :: ib

    normint_mass = (mrn_amax**(mrn_q+4.0D0) - mrn_amin**(mrn_q+4.0D0)) / (mrn_q+4.0D0)
    A_mrn = (1.0D-2*metallicity)*density*MH / ((4.0D0/3.0D0)*PI*mrn_rhogr*normint_mass)

    G = 0.0D0
    DO ib=1,mrn_nbin
       bin_a  = mrn_amin*(mrn_amax/mrn_amin)**((REAL(ib,KIND=dp)-0.5D0)/REAL(mrn_nbin,KIND=dp))
       bin_da = bin_a*LOG(mrn_amax/mrn_amin)/REAL(mrn_nbin,KIND=dp)
       tau    = bin_a*KB*temperature/(REAL(zion,KIND=dp)**2*EC*EC)
       Jt     = 1.0D0 + SQRT(PI/(2.0D0*tau))
       G = G + A_mrn*bin_a**mrn_q * PI*bin_a**2 * Jt * bin_da
    ENDDO

  END FUNCTION MRN_DS_G

END MODULE grain_distribution_module
#endif
