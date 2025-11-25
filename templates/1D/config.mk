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
RESTART              = 0
OUTRAYINFO           = 0
CHEMANALYSIS         = 0
