CMP     =  gcc
DEFS_ADD = 

exeName =  ColumeDensity
dir = 
inc = 
lib = 

ifeq ($(CMP),intel)
  FortC = mpiifort
  CFLAG = -fpp -O3 -mcmodel=large -fp-model fast=2 -g#-heap-arrays 20#-ipo -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2
else ifeq ($(CMP),gcc)
  FortC = gfortran
  CFLAG = -cpp -O3 -Warray-bounds -fbacktrace -g -ffree-line-length-none -mcmodel=large -funroll-loops -floop-optimize
endif

SrcT   := m_Mesh.f90 m_Ray_box.f90 init3DPDR.f90
src:= $(addprefix $(dir), ${SrcT})

all: $(exeName)
src_obj  = $(src:%.f90=%.o)
$(exeName):$(src_obj)
	$(FortC) $(CFLAG) -fopenmp -o $@ $(src_obj) $(lib)
$(src_obj):%.o :%.f90
	$(FortC) $(CFLAG) -fopenmp $(inc) $(DEFS_ADD) -c $<

.PHONY: clean
clean:
	rm -fr *.o *.mod $(exeName) $(dir)*.o

