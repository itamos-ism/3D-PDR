subroutine allocations
use healpix_types
use maincode_module
use global_module
use uclpdr_module
!use maincode_local, only : prev_cooling

!load SCO_GRID data 
SCO_GRID(1:8,1) = (/0.000D+00,-1.408D-02,-1.099D-01,-4.400D-01,-1.154D+00,-1.888D+00,-2.760D+00,-4.001D+00/)
SCO_GRID(1:8,2) = (/-8.539D-02,-1.015D-01,-2.104D-01,-5.608D-01,-1.272D+00,-1.973D+00,-2.818D+00,-4.055D+00/)
SCO_GRID(1:8,3) = (/-1.451D-01,-1.612D-01,-2.708D-01,-6.273D-01,-1.355D+00,-2.057D+00,-2.902D+00,-4.122D+00/)
SCO_GRID(1:8,4) = (/-4.559D-01,-4.666D-01,-5.432D-01,-8.665D-01,-1.602D+00,-2.303D+00,-3.146D+00,-4.421D+00/)
SCO_GRID(1:8,5) = (/-1.303D+00,-1.312D+00,-1.367D+00,-1.676D+00,-2.305D+00,-3.034D+00,-3.758D+00,-5.077D+00/)
SCO_GRID(1:8,6) = (/-3.883D+00,-3.888D+00,-3.936D+00,-4.197D+00,-4.739D+00,-5.165D+00,-5.441D+00,-6.446D+00/)


allocate(species(1:nspec))
allocate(init_abundance(1:nspec))
allocate(mass(1:nspec))
allocate(reactant(1:nreac,1:3))
allocate(product(1:nreac,1:4))
allocate(rate(1:nreac))
allocate(alpha(1:nreac))
allocate(beta(1:nreac))
allocate(gamma(1:nreac))
allocate(rtmin(1:nreac))
allocate(rtmax(1:nreac))
allocate(duplicate(1:nreac))
#ifdef CHEMANALYSIS
allocate(temp_rate(1:nreac,1:pdr_ptot))
#endif


return
end subroutine allocations


