---
title: '3D-PDR (2025): modelling three-dimensional photodissociation regions'
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
  - name: Xuefei Tang
    orcid:
    equal-contrib: true
    affiliation: 1
  - name: Xue-Jian Jiang
    orcid: 0000-0002-8899-4673
    equal-contrib: true
    affiliation: 1
  - name: Gaoyuan Zhang
    orcid: 
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Research Center for Astronomical Computing, Zhejiang Lab, 311100, Hangzhou, China
   index: 1
 - name: Faculty of Physics, University of Duisburg-Essen, Lotharstraß 1, 47057, Duisburg, Germany
   index: 2
date: June 2026
bibliography: 3DPDR_joss.bib
---

# Summary

Photodissociation regions (PDR) are dynamic interfaces in the interstellar medium (ISM) where ultraviolet photons (UV) from nearby massive stars strongly influence the physical and chemical state of the gas and dust [@Tielens1985; @Hollenbach1999; @Wolfire2022]. PDRs are regions acting as transition zones between the ionized and molecular phases and where the atomic-to-molecular phase is located. They are characterized by the dissociation of molecules, especially that of $\rm H_2$, as well as the ionization of atoms, leading to the formation of complex structures and the emission of molecular, atomic and ionic lines. PDRs are found in environments ranging from molecular cloud edges to surfaces of star-forming regions, as well as HII regions. Understanding PDRs is crucial for studying the life-cycle of the ISM, the environmental conditions leading to star-formation, and consequently the evolution of galaxies.

Numerical modeling of PDRs presents significant challenges, due to the coupling of the given model distribution with radiative transfer and chemical kinetics involving hundreds to thousands of reactions in non-equilibrium environments. Various models of PDRs are often one-dimensional, assuming simple slab or spherical geometries [@Roellig2007]. While computationally efficient, such models cannot capture the inherently complex, three-dimensional structures of real astrophysical environments that shape the chemistry and emission properties. A fully three-dimensional approach is, therefore, essential to realistically model these environments, allowing researchers to connect theory with the detailed observations now possible with radio-telescopes.

The 3D-PDR code has been one of the earliest attempts in modelling the PDR chemistry in three-dimensional density distributions [@Bisbas2012]. The upgraded 3D-PDR, offers a significant speed-up in calculations of three-dimensional structures, thanks to the new ray-tracing algorithm RAYTHEIA [@Zhu2026]. It performs thermal balance iterations and terminates once the total heating matches with the total cooling to within a user-defined tolerance parameter. It outputs the abundance distribution of species of the given chemical network, gas and dust temperatures, local heating and cooling for all individual processes, as well as emissivities and level populations of the coolants considered. The latter can be further post-processes to construct synthetic observations by solving the radiative transfer equation from the point of view of the observer. 

# Statement of need

The 3D-PDR code provides a state-of-the-art numerical framework to simulate the chemistry, thermal balance, and line emission of PDRs. It can be used as a one-dimensional code for fast PDR calculations [@Bisbas2014; @Bisbas2015; @Gaches2015; @Gaches2019a; @Gaches2019b; @Luo2023; @Luo2024] and the construction of thermochemical grids [@Bisbas2019; @Bisbas2023; @Bisbas2025; @Dasyra2022], and also as a three-dimensional tool to post-process snapshots of molecular clouds and star-forming regions modelled with hydrodynamical codes [@Bisbas2017a; @Bisbas2017b; @Bisbas2018; @Bisbas2021; @Bisbas2024; @Gaches2022a; @Gaches2022b; @Obolentseva2024]. To our knowledge, 3D-PDR is the only available code for treating PDRs in three-dimensions. The upgraded version of 3D-PDR presented here includes significant improvements in computational performance, modularity, and usability, enabling large-scale parameter studies and integration with hydrodynamical processes. In doing so, it advances the field beyond earlier one-dimensional codes and provides an important tool for researchers seeking to interpret observations in the era of facilities such as ALMA and JWST.

# Software design 

The upgraded 3D-PDR code features the new ray-tracing algorithm RAYTHEIA [@Zhu2026]. RAYTHEIA is a high-performance reverse ray-tracing algorithm designed to efficiently solve three-dimensional direction-dependent equations in astrophysical simulations, with particular focus on photodissociation regions (PDRs). The method adopts a dual-grid framework in which the native simulation mesh serves as the source grid for ray emission, while an adaptive Cartesian contribution grid is used for efficient ray traversal and accumulation of directional quantities. Rays are discretized using the HEALPix scheme and treated as infinitesimally thin pencil beams.

Key innovations of RAYTHEIA include a leaf-only linear-octree data structure combined with Morton-code indexing and digital differential analyser (DDA) traversal for efficient ray-walking, an analytical slab method for exact path-length calculations through grid cells, and a hybrid MPI/OpenMP distributed-memory framework employing a chunk-to-chunk communication strategy. These advances provide near-ideal parallel scalability, substantially reduce memory requirements, and enable high-resolution three-dimensional PDR simulations. RAYTHEIA is modular, computationally efficient, and can be readily integrated into existing simulation codes requiring fast solutions of direction-dependent radiative-transfer problems.

# Research impact statement

The upgraded 3D-PDR code, powered by the new RAYTHEIA ray-tracing algorithm, enables three-dimensional astrochemical and radiative-transfer calculations at resolutions that were previously computationally prohibitive. By combining a dual-grid framework, adaptive octree traversal, analytical path-length calculations, and highly scalable hybrid parallelization, the code can model the chemistry, thermal balance, and line emission of PDRs in high-resolution simulations. This capability bridges the gap between modern large-scale hydrodynamical simulations and detailed astrochemical modelling, allowing researchers to study the interaction between radiation and the interstellar medium in three-dimensional distributions.

The software provides a versatile platform for investigating a wide range of astrochemical problems, including diffuse and molecular cloud evolution in different ISM environments, star-forming regions, CO-dark molecular gas, and the interpretation of observations from archival data, flagship facilities such as ALMA and JWST, and future instruments. Through its modular design, RAYTHEIA can be integrated into existing numerical frameworks requiring efficient solutions of direction-dependent equations, extending its applicability beyond PDR modelling. The combination of high performance and scalability makes the upgraded 3D-PDR framework a valuable tool for both theoretical studies and the generation of synthetic observations used to connect simulations with observed data cubes.


# Example case

As an example, we post-process a sub-region from the SILCC-Zoom hydrodynamical simulation of @Seifried2017 modelling a star-forming molecular cloud. The selected region has a mass of $M_{\rm tot}=7.3\times10^4\,{\rm M}_{\odot}$ and an extent of 62.5 pc resolved with $256^3$ uniform cells achieving a spatial resolution of 0.24 pc. It contains a filamentary structure with significant substructure, including dense clumps and extended low-density regions. We consider an FUV strength of $\chi/\chi_0=1$ [@Draine1978] and the $\cal H$-function of @Padovani2018 describing the cosmic-ray ionization rate and its attenuation along depth. We assume solar metallicity, dust-to-gas ratio of $10^{-2}$ and a sub-set of the UMIST2012 chemical network [@McElroy2013] consisting of 77 species and 1158 reactions. 

Figure \@ref(fig:example) shows example outputs of the velocity integrated emission maps of \[CII\] 158 μm, \[CI\] (1-0), CO $J=1-0$, \[OI\] 63 μm, HCO$^+$ (1-0) and HCN (1-0). Similarly, the column densities of all modelled species can be plotted. 

![Velocity integrated emission maps. Top row shows \[CII\] 158 μm, \[CI\](1-0) and CO(1-0). Bottom row shows \[OI\] 63 μm, HCO+ (1-0) and HCN (1-0) emission lines.](six_line_maps_2x3.png){#fig:example width=100%}.


# Acknowledgements

This work is supported by the Leading Innovation and Enterpreneurship Team of Zhejiang Province of China (Grant No. 2023R01008). BALG is supported by the German Research Foundation (DFG) in the form of an Emmy Noether Research Group - DFG project
#542802847 (GA 3170/3-1).

# References
