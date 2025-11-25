include config.mk

#SUNDIALS flags
INCLUDES          = -I${SUNDIALS_DIR}/include
LIBRARIES	  = -L${SUNDIALS_DIR}/lib
LIBS	          = -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_sunnonlinsolnewton -lsundials_core -lm


ifneq ($(DIMENSIONS),1)
  ifneq ($(DIMENSIONS),3)
    $(error "Invalid value for DIMENSIONS. Only 1 and 3 are allowed")
  endif
endif

ifeq ($(DIMENSIONS),1)
  ifeq ($(RAYTHEIA),1)
      $(error "Cannot compile with DIMENSIONS = 1 and RAYTHEIA = 1")
  endif
  ifeq ($(XYZ),1)
      $(error "Cannot compile with DIMENSIONS = 1 and XYZ = 1")
  endif
endif

ifeq ($(THERMALBALANCE),0)
  ifeq ($(GUESS_TEMP),1)
      $(error "Cannot compile with THERMALBALANCE = 0 and GUESS_TEMP = 1")
  endif
endif

ifeq ($(PYWRAP),1)
  ifeq ($(DIMENSIONS),3)
      $(error "PYWRAP output not compatible with DIMENSIONS = 3")
  endif
endif

#F95 compiler
ifeq ($(F90),f95) 
  ifeq ($(OPENMP),1)
    OPT += -fopenmp -DOPENMP
    OPTC += -fopenmp
  endif
#ifort compiler
else ifeq ($(F90),ifort)
  ifeq ($(OPENMP),1) 
    CFLAGS += -DOPENMP
    OPT += -openmp
    OPTC += -fopenmp
  endif
#gfortran compiler
else ifeq ($(F90),gfortran)
  ifeq ($(OPENMP),1)
    CFLAGS += -DOPENMP
    OPT += -fopenmp
    OPTC += -fopenmp
  endif
#gfortran44 compiler
else ifeq ($(F90),gfortran44)
  ifeq ($(OPENMP),1)
    CFLAGS += -DOPENMP
    OPT += -fopenmp
    OPTC += -fopenmp
  endif
endif
ifeq ($(OPTIMISE),0)
  OPT += -O0 -fcheck=all -Wall
  OPTC += -O0 -Wall
endif
ifeq ($(OPTIMISE),1)
  OPT += -O1
  OPTC += -O1
endif
ifeq ($(OPTIMISE),2)
  OPT += -O2 
  OPTC += -O2
endif
ifeq ($(OPTIMISE),3)
  OPT += -O3
  OPTC += -O3
endif
ifeq ($(PYWRAP),1)
  CFLAGS += -DPYWRAP
endif

#Chemistry flags
ifeq ($(THERMALBALANCE),1)
  CFLAGS += -DTHERMALBALANCE
endif
ifeq ($(GUESS_TEMP),1)
  CFLAGS += -DGUESS_TEMP
endif
ifeq ($(DIMENSIONS),1)
  CFLAGS += -DONEDIMENSIONAL
endif

ifeq ($(DUST),HTT91)
  CFLAGS += -DDUST
endif
ifeq ($(FORCECONVERGENCE),1)
  CFLAGS += -DFORCECONVERGENCE
endif
ifeq ($(GRAINRECOMB),1)
  CFLAGS += -DGRAINRECOMB
endif
ifeq ($(SUPRATHERMAL),1)
  CFLAGS += -DSUPRATHERMAL
endif
ifeq ($(H2FORM),CT02)
  CFLAGS += -DH2FORM_CT02
endif
ifeq ($(H2FORM),R07)
  CFLAGS += -DH2FORM_R07
endif
ifeq ($(H2FORM),SIMPLE)
  CFLAGS += -DH2FORM_SIMPLE
endif

ifeq ($(NETWORK),REDUCED)
  CFLAGS += -DREDUCED
endif
ifeq ($(NETWORK),MEDIUM)
  CFLAGS += -DMEDIUM
endif
ifeq ($(NETWORK),FULL)
  CFLAGS += -DFULL
endif
ifeq ($(NETWORK),MYNETWORK)
  CFLAGS += -DMYNETWORK
endif
ifeq ($(RAYTHEIA),1)
  CFLAGS += -DRAYTHEIA
  CFLAGS += -DRAYTHEIA_MO
endif
ifeq ($(RAYTHEIA),2)
  CFLAGS += -DRAYTHEIA
endif
ifeq ($(XYZ),1)
  CFLAGS += -DXYZ
endif
ifneq ($(CRATTENUATION),0)
  CFLAGS += -DCRATTENUATION=$(CRATTENUATION)
endif
ifeq ($(RESTART),1)
  CFLAGS += -DRESTART
endif
ifeq ($(OUTRAYINFO),1)
  CFLAGS += -DOUTRAYINFO
endif
ifeq ($(CHEMANALYSIS),1)
  CFLAGS += -DCHEMANALYSIS
endif


MODULE_OBJ += definitions.o healpix_types.o modules.o m_Mesh.o
CODE_OBJ += healpix.o input_parameters.o solvlevpop.o \
read_species.o read_rates.o heapsort.o calc_reac_rates.o shield.o \
spline.o m_Ray_box.o escape_probability.o eval_points.o heatingfunctions.o \
find_Ccoeff.o read_input.o read_command_line.o \
logo.o writeparams.o readdensity.o read_coolants.o initialization.o \
changetemperature.o UVfield.o CRattenuation.o columndensity.o partition_function.o \
lte_population.o allocations.o writeoutputs.o chemicaliterations.o \
coolingfunctions.o checkconvergence.o
CODE_OBJ += calculate_abundances.o
ifeq ($(CHEMANALYSIS),1)
CODE_OBJ += analyse_chem.o
endif



ifeq ($(NETWORK),REDUCED)
CODE_OBJ += odes_reduced.o
endif
ifeq ($(NETWORK),MEDIUM)
CODE_OBJ += odes_medium.o
endif
ifeq ($(NETWORK),FULL)
CODE_OBJ += odes_full.o
endif
ifeq ($(NETWORK),MYNETWORK)
CODE_OBJ += odes_mynetwork.o
endif
ifeq ($(DUST),HTT91)
CODE_OBJ += dust_t.o
endif
ifeq ($(H2FORM),CT02)
CODE_OBJ += h2_form.o
endif

OBJ += $(MODULE_OBJ) $(CODE_OBJ)
#OPT += -fcheck=all

%.o: %.F90 definitions.o healpix_types.o modules.o m_Mesh.o
	$(F90) $(CPPFLAGS) $(OPT) $(CFLAGS) -c $<

%.o: %.F
	$(F90) $(OPT) $(CFLAGS) -c $<

%.o: %.c
	$(CC) $(OPTC) $(INCLUDES) $(CFLAGS) -c $<

3DPDR :: $(OBJ) 3DPDR.o
	$(F90) $(OPT) $(CFLAGS) $(LIBRARIES) -o 3DPDR $(OBJ) 3DPDR.o $(LIBS)
	mv 3DPDR ../
clean :: 
	rm *.mod *.o ../restart.bin

definitions.o : definitions.F90
	$(F90) $(CPPFLAGS) $(OPT) $(CFLAGS) -c $<

healpix_types.o : healpix_types.F90
	$(F90) $(CPPFLAGS) $(OPT) $(CFLAGS) -c $<

modules.o : modules.F90
	$(F90) $(CPPFLAGS) $(OPT) $(CFLAGS) -c $<

m_Mesh.o : m_Mesh.F90
	$(F90) $(CPPFLAGS) $(OPT) $(CFLAGS) -c $<

healpix_types.o : definitions.o
modules.o : definitions.o healpix_types.o
