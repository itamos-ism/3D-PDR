F90 = ifx
HDF5_DIR = /mnt/d/linux_workspace/software/hdf5
HDF5_INCLUDE = $(HDF5_DIR)/include
HDF5_LIB = $(HDF5_DIR)/lib
CPPFLAGS += -I$(HDF5_INCLUDE)
LIBRARIES += -L$(HDF5_LIB)
LIBS     = -lhdf5_fortran -lhdf5

OPT = -O3
OBJ = m_paramters.o m_visual.o m_writeoutputs.o

%.o: %.F90
	$(F90) $(CPPFLAGS) $(OPT) -c $<

convert :: $(OBJ) convert_fin2h5.o
	$(F90) $(OPT) $(LIBRARIES) -o convert $(OBJ) convert_fin2h5.o $(LIBS)
clean :: 
	rm convert *.mod *.o *.xmf
