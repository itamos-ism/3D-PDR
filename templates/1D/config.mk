#==================================================================
#   	 	       3D-PDR Makefile
#----------------Compiler options----------------------------------
#OPENMP   	  : 1 - run in parallel (OpenMP)
#         	    0 - run in serial (1 CPU) 
#OPTIMISE 	  : use 0 or 1 for development & debugging
#     	            use 2 or 3 for running models
#PYWRAP		  : 1 - output for python wrapper
#		  : 0 - switch off
#-----------------Chemistry options--------------------------------
#DIMENSIONS 	  : 1 - compile for 1D models (needs RAYTHEIA = 0)
#		    3 - compile for 3D models
#RAYTHEIA 	  : 1 - switch on with memory optimization (slower code, less RAM, recommended)
#	 	  : 2 - switch on without memory optimization (faster code, more RAM)
#	 	  : 0 - switch off (uses the 2012 method, slower code, high RAM, no grid necessary)
#XYZ              : 0 - z coordinate changes first - xyz
#	   	  : 1 - x coordinate changes first - xyz
#NETWORK 	  : REDUCED (33 species & 331 reactions)
#         	    MEDIUM  (77 species & 1158 reactions; due to BALG)
#         	    FULL    (215 species & 2926 rections)
#DUST   	  : HTT91 - Hollenbach, Takahashi & Tielens (1991)
#                   0     - for isothermal dust models
#THERMALBALANCE   : 1 - for thermal balance iterations
#		    0 - for isothermal runs (needs GUESS_TEMP=0 also)
#FORCECONVERGENCE : 1 - Helps fast convergence (recommended)
#		    0 - Switch off
#GRAINRECOMB      : 1 - Electron recombination on dust grains
#	            0 - Switch off
#SUPRATHERMAL     : 1 - Suprathermal formation of CO via CH+
#	            0 - Switch off
#H2FORM	          : CT02   - Cazaux & Tielens (2002,2004) treatment
#	            SIMPLE - 3e-18*sqrt(Tgas)*exp(-Tgas/1e3) expression
#	            R07    - Roellig+07 benchmarking paper
#CRATTENUATION    : 1 - models CR attenuation (L/H of Padovani+ models)
#                   0 - switch off (default)
#RESTART          : 1 - restarts an interrupted model
#                   0 - switch off (default)
#OUTRAYINFO       : 1 - saves info for each ray in different files
#		  : 0 - switch off (default)
#CHEMANALYSIS     : 1 - outputs chemical analysis per point (not recommended for 3D models)
#		  : 0 - switch off (default)
#CHEMSTEADY       : 1 - Continue each chemical integration in doubling time
#                       blocks (up to 128x the evolution time) until no
#                       species changes by >1% over a doubling. Otherwise the
#                       chemical age of a grid point depends on how many
#                       chemistry calls it received before its thermal balance
#                       converged, leaving unphysical abundance jumps between
#                       neighbouring points (recommended)
#                   0 - Switch off (each call advances by the evolution
#                       time only)
#CHEMRETRY        : 1 - Re-integrate a grid cell whose CVODE chemistry solve
#                       fails (e.g. "error test failed ... |h|=hmin at t=0"
#                       for very stiff chemistry under extreme conditions) from
#                       the saved initial abundances with a forced, shrinking
#                       initial step (up to 4 retries). Recovers the few stiff
#                       cells that CVODE's automatic initial step over-reaches;
#                       the first attempt is unchanged so normal cells are
#                       unaffected (recommended)
#                   0 - Switch off (failed cells abandoned at best effort)
#NGACCEL          : 1 - Ng (1974) acceleration of the level population
#                       iterations (faster convergence, recommended)
#                   0 - Switch off
#NGCYCLE          : 1 - Dynamic per-point damping of oscillating (flip-flop)
#                       level populations: whenever an update reverses the
#                       previous one, the last two iterates are averaged.
#                       Acts on each point/coolant individually, every
#                       iteration, as soon as the oscillation appears,
#                       instead of stalling until the FORCECONVERGENCE
#                       averaging at level iteration 75
#                       (needs NGACCEL=1; recommended)
#NGRELAX          : 1 - Adaptive per-point under-relaxation of persistently
#                       oscillating level populations: a per-cell counter
#                       shrinks the relaxation weight (1/(2+n)) while a cell
#                       keeps reversing, damping the violent CO(1-0) flip-flop
#                       near the inversion boundary (excitation temp -> inf)
#                       that fixed weight-1/2 cannot, in hot CR-heated cells
#                       at n~n_crit. Needs NGACCEL=1; supersedes the NGCYCLE
#                       damping when both are on (recommended, esp. for 3D)
#                   0 - Switch off
#ILLINOIS         : 1 - Illinois-modified regula-falsi (in ln(T)) for the
#                       thermal-balance temperature refinement, replacing
#                       plain bisection (faster convergence, recommended)
#                   0 - Switch off (use plain bisection)
#REBRACKET        : 1 - When the thermal-balance search stalls with a large
#                       residual imbalance (bracket confined away from the
#                       root by early, not-yet-relaxed cooling rates),
#                       re-open the bracket and resume the search (up to 3
#                       attempts per point) instead of forcing convergence
#                       (recommended)
#                   0 - Switch off (force convergence on stall)
#SOBOLEV          : 1 - Sobolev (ALI-like) net radiative rates: A*beta and
#                       beta*B*BB instead of the mean field (1-beta)*S+beta*BB.
#                       Same solution, but removes the lagged self-coupling
#                       that makes optically thick lines falsely converge
#                       (<1% change per iteration while far from equilibrium),
#                       which leaves spurious cooling at low Tgas (recommended)
#                   0 - Switch off (original mean-field iteration)
#------------------------------------------------------------------
F90                  = gfortran
CC                   = gcc
CPPFLAGS             = -cpp
OPENMP               = 1
OPTIMISE             = 3
PYWRAP               = 0
DIMENSIONS           = 1
RAYTHEIA             = 0
XYZ                  = 0
NETWORK              = REDUCED
DUST                 = HTT91
GUESS_TEMP           = 1
THERMALBALANCE       = 1
FORCECONVERGENCE     = 1
GRAINRECOMB          = 0
SUPRATHERMAL         = 0
H2FORM               = CT02
CRATTENUATION        = 0
#Convergence of ODEs
CHEMSTEADY           = 1
CHEMRETRY            = 1
#Ng acceleration
NGACCEL              = 1
NGCYCLE              = 1
NGRELAX              = 1
#Thermal balance
ILLINOIS             = 1
REBRACKET            = 1
#Escape probability
SOBOLEV              = 1
###
RESTART              = 0
OUTRAYINFO           = 0
CHEMANALYSIS         = 0
