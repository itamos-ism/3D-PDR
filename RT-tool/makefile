#--------RT tool makefile-------
#--------Thomas G. Bisbas-------

FC = gfortran
FFLAGS = -g -fbacktrace -fcheck=all
DUST = 0
#LEVMAX : 1 - read the LMAX= header line written by 3D-PDR when compiled with LEVMAX=1
#             and truncate each coolant to the same level cap used in the 3D-PDR run.
#             Must match the 3D-PDR build (both LEVMAX=1, or both LEVMAX=0).
#         0 - switch off (default; .RTspop.fin files from LEVMAX=0 runs)
LEVMAX = 0

ifeq ($(DUST),1)
   OPT += -DDUST
endif
ifeq ($(LEVMAX),1)
   OPT += -DLEVMAX
endif

SRC = modules.F90 readcoolants.F90 read_input.F90 readparams.F90 \
      readpdrfiles.F90 initializations.F90 writeoutputs.F90 excitation.F90 \
      readspopfile.F90 opticaldepth.F90 radiativetransfer.F90 RTtool.F90

EXEC = RTtool

all: $(EXEC)
	mv $(EXEC) ../

$(EXEC): $(SRC)
	$(FC) $(FFLAGS) $(OPT) -o $(EXEC) $(SRC)

clean:
	rm -f $(EXEC) *.o *.mod

