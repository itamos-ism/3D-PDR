---
title: 'Modelling photodissociation regions with the 3D-PDR code'
authors:
  - name: Thomas G. Bisbas
    orcid: 0000-0003-2733-4580
    equal-contrib: true
    affiliation: 1 # ("1,2" Multiple affiliations must be quoted)
  - name: Zhengping Zhu
    orcid:
    equal-contrib: true
    affiliation: 1
  - name: Brandt A. L. Gaches
    orcid: 0000-0003-4224-6829
    equal-contrib: true
    affiliation: 2
  - name: Xue-Jian Jiang
    orcid: 0000-0002-8899-4673
    equal-contrib: true
    affiliation: 1
  - name: Xuefei Tang
    orcid:
    equal-contrib: true
    affiliation: 1
  - name: Gaoyuan Zhang
    orcid: 
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Research Center for Computational Earth and Space Science, Zhejiang Lab, 311100, Hangzhou, China
   index: 1
 - name: Faculty of Physics, University of Duisburg-Essen, Lotharstraß 1, 47057, Duisburg, Germany
   index: 2
date: June 2026
bibliography: paper.bib
---

# Summary

3D-PDR is an open-source astrochemical code designed to model the chemistry, thermal balance, and line emission of photodissociation regions (PDRs) in one-, and three-dimensional environments. The code solves the coupled interaction between far-ultraviolet, gas-phase chemistry, dust processes, thermal physics, and level populations, enabling the prediction of molecular abundances, temperatures, and observable emission from the interstellar medium. The upgraded version presented here introduces the high-performance ray-tracing algorithm RAYTHEIA, together with new physical modules and a modular software architecture that significantly improve computational efficiency, scalability, and scientific capabilities. The code provides predictions using 1D and 3D density distributions, including post-processing of hydrodynamical simulations and the generation of synthetic observations for comparison with data from current and future instruments.


# State of the field

Photodissociation regions are dynamic interfaces in the interstellar medium (ISM) where ultraviolet photons (UV) from nearby massive stars strongly influence the physical and chemical state of the gas and dust [@Tielens1985; @Hollenbach1999; @Wolfire2022]. PDRs are regions acting as transition zones between the ionized and molecular phases and where the atomic-to-molecular phase is located. They are characterized by the dissociation of molecules, especially that of $\rm H_2$, as well as the ionization of atoms, leading to the formation of complex structures and the emission of molecular, atomic and ionic lines. PDRs are found in environments ranging from molecular cloud edges to surfaces of star-forming regions, as well as HII regions. Since the bulk of the neutral ISM emission is associated with PDRs [@Hollenbach1999; @Wolfire2022], understanding them is crucial for studying the small and large scale of the ISM, including the evolution of galaxies.

Numerical modeling of PDRs presents significant challenges, due to the coupling of the given model distribution with radiative transfer and chemical pathways involving hundreds to thousands of reactions in non-equilibrium environments. Various models of PDRs are often one-dimensional, assuming simple slab or spherical geometries [@Roellig2007]. 

The 3D-PDR code has been one of the earliest attempts in modelling the PDR chemistry in three-dimensional density distributions [@Bisbas2012]. In particular, 3D-PDR can be used as a one-dimensional code for fast PDR calculations [@Bisbas2014; @Bisbas2015; @Gaches2019a; @Gaches2019b; @Luo2023; @Luo2024] and the construction of one-dimensional thermochemical grids covering a wide range of ISM environmental parameters [@Gaches2015; @Bisbas2019; @Bisbas2023; @Bisbas2025; @Dasyra2022]. It can also be used as a three-dimensional tool to post-process snapshots of molecular clouds and star-forming regions modelled with hydrodynamical codes [@Bisbas2017a; @Bisbas2017b; @Bisbas2018; @Bisbas2021; @Bisbas2024; @Gaches2022a; @Gaches2022b; @Obolentseva2024] and to generate synthetic observations for direct comparison with observational data.

The upgraded 3D-PDR presented here, offers a significant speed-up in calculations of three-dimensional structures, thanks to the new ray-tracing algorithm RAYTHEIA [@Zhu2026]. 

# Statement of need

While computationally efficient with highly sophisticated chemistry, one-dimensional models cannot capture the complexity of three-dimensional structures that real interstellar clouds exhibit and therefore their abundances distribution and emission properties. A fully three-dimensional approach is, therefore, essential to realistically model these environments, allowing researchers to connect theory with the detailed observations now possible with radio-telescopes.

The 3D-PDR code provides a state-of-the-art numerical framework to simulate the chemistry, thermal balance, and line emission of PDRs. To our knowledge, 3D-PDR is the only available code for treating PDRs in three-dimensions. The upgraded version of 3D-PDR includes significant improvements in computational performance, modularity, and usability, enabling large-scale parameter studies and integration with hydrodynamical processes. In doing so, it advances the field beyond earlier one-dimensional codes and provides an important tool for researchers seeking to interpret observations in the era of facilities such as ALMA and JWST.

# Research impact statement

The upgraded 3D-PDR code enables three-dimensional astrochemical and radiative-transfer calculations at resolutions that were previously computationally prohibitive. By combining a dual-grid framework, adaptive octree traversal, analytical path-length calculations, and efficient hybrid parallelization, the code can model the chemistry, thermal balance, and line emission of PDRs in high-resolution three-dimensional simulations. This capability bridges the gap between modern large-scale hydrodynamical simulations and detailed astrochemical modelling, allowing the community to study the interaction between radiation and the ISM in full three-dimensions.

The software provides a platform for investigating a broad range of astrochemical problems, including molecular cloud evolution, star formation, feedback processes, and the interpretation of observations from flagship facilities such as ALMA and JWST, and future instruments. This makes 3D-PDR a powerful tool for both theoretical studies and the generation of synthetic observations used to connect simulations with observed data.


# Software design

The upgraded 3D-PDR code features the new ray-tracing algorithm RAYTHEIA [@Zhu2026]. RAYTHEIA is a high-performance reverse ray-tracing algorithm designed to efficiently solve three-dimensional direction-dependent equations in astrophysical simulations, with particular focus on photodissociation regions (PDRs). The method adopts a dual-grid framework in which the native simulation mesh serves as the source grid for ray emission, while an adaptive Cartesian contribution grid is used for efficient ray traversal and accumulation of directional quantities. Rays are discretized using the HEALPix scheme [@Gorski2005] and treated as infinitesimally thin pencil beams.

Key innovations of RAYTHEIA include a leaf-only linear-octree data structure combined with Morton-code indexing and digital differential analyser (DDA) traversal for efficient ray-walking, an analytical slab method for exact path-length calculations through grid cells, and a hybrid MPI/OpenMP distributed-memory framework employing a chunk-to-chunk communication strategy. These advances provide near-ideal parallel scalability, substantially reduce memory requirements, and enable high-resolution three-dimensional PDR simulations. RAYTHEIA is modular, computationally efficient, and can be readily integrated into existing simulation codes requiring fast solutions of direction-dependent radiative-transfer problems.

3D-PDR performs thermal balance and level populations iterations and terminates once the total heating matches with the total cooling to within a user-defined tolerance parameter. It outputs the abundance distribution of species of the given chemical network, gas and dust temperatures, local heating and cooling for all individual processes, as well as emissivities and level populations of the coolants considered. The latter can be further post-processes to construct synthetic observations by solving the radiative transfer equation from the point of view of the observer. 

The code includes several new physical modules from those initially presented in @Bisbas2012 that extend its applicability to a wider range of interstellar environments. These additions include routines for modelling cosmic-ray energy spectra and their attenuation through gas columns [@Gaches2019a; @Gaches2022a], allowing for a more accurate treatment of cosmic-ray ionization rates. The code also incorporates the suprathermal formation pathway of CO via $\mathrm{CH}^+$ [@Visser2009; @Bisbas2019], which becomes important in regions where non-thermal processes enhance the formation of CO at low column densities. Furthermore, the chemical network has been expanded to include electron-ion recombination on dust grain surfaces [@Weingartner2001], providing a better description of charge balance in dense and shielded gas. Finally, the upgraded version introduces the calculation of H$_2$ ro-vibrational and rotational line emission (Gaches et al., submitted), enabling direct comparison with JWST data. 


# Example case

As an example, we post-process a sub-region from the SILCC-Zoom hydrodynamical simulation of @Seifried2017 modelling a star-forming molecular cloud. The selected region has a mass of $M_{\rm tot}=7.3\times10^4\,{\rm M}_{\odot}$ and an extent of 62.5 pc resolved with $256^3$ uniform cells achieving a spatial resolution of 0.24 pc. It contains a filamentary structure with significant substructure, including dense clumps and extended low-density regions. We consider an FUV strength of $\chi/\chi_0=10$ [@Draine1978] and the $\cal L$-function of @Padovani2018 describing the cosmic-ray ionization rate and its attenuation along depth. We assume solar metallicity, dust-to-gas ratio of $10^{-2}$ and a sub-set of the UMIST2012 chemical network [@McElroy2013] consisting of 79 species and 1169 reactions. The initial elemental abundances were set as follows: $\mathrm{C}^+=10^{-4}$, $\mathrm{O}=3\times10^{-4}$, $\mathrm{N}=3.3\times10^{-5}$, $\mathrm{Mg}^+=2.7\times10^{-7}$, and $\mathrm{He}=8.5\times10^{-2}$. 

![Top half: column densities of $\mathrm{C}^+$, C and CO (upper row) and O, $\mathrm{HCO}^+$ and HCN (lower row). Bottom half: the corresponding velocity-integrated emission maps, [CII] 158 μm, [CI] (1-0) and CO (1-0) (upper row), and [OI] 63 μm, $\mathrm{HCO}^+$ (1-0) and HCN (1-0) (lower row). Note that the latter two lines are in log-scale. \label{fig:example}](six_emission_and_cds.png){ width=95% }

\autoref{fig:example} shows example outputs of the column densities and velocity integrated emission maps. The top half shows column densities of $\mathrm{C}^+$, C and CO (upper row) and O, $\mathrm{HCO}^+$ and HCN (lower row). The bottom half shows the corresponding velocity-integrated emission maps, [CII] 158 μm, [CI] (1-0) and CO (1-0) (upper row), and [OI] 63 μm, $\mathrm{HCO}^+$ (1-0) and HCN (1-0) (lower row). Note that the colour bars of the latter two velocity integrated emission maps are in log-scale. 



<!-- 
Figure \@ref(fig:example) shows example outputs including the column densities of ${\rm H}_2$, $\rm C^+$, C, and CO, as well as the velocity integrated emission maps of \[CII\] 158 μm, \[CI\] (1-0) and CO $J=1-0$. 
![The left panel shows the H$_2$ column density. Top row shows the velocity integrated maps of \[CII\] 158$\mu$m, \[CI\](1-0) and CO(1-0). Bottom row shows the corresponding column densities.](figure_joss.png){#fig:example width=100%}.
-->

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing of this manuscript, or the preparation of supporting materials.

# Acknowledgements

This work is supported by the Leading Innovation and Enterpreneurship Team of Zhejiang Province of China (Grant No. 2023R01008). BALG is supported by the German Research Foundation (DFG) in the form of an Emmy Noether Research Group - DFG project #542802847 (GA 3170/3-1).

# References
