# ifndef FOLDER
# 	FOLDER="../"
# endif
#random_folder=${FOLDER}/fortran_rng

# Makefile

FC = ifx
FFLAGS= -i8 -warn all -Ofast -static
# -g -fcheck=all -fbacktrace
STATICLIB= /usr/local/lib/libarpackILP64.a  ${MKLROOT}/lib/libmkl_blas95_ilp64.a ${MKLROOT}/lib/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
#MKL_LINK=-L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
MKL_INCLUDE=-I"${MKLROOT}/include"
EXE_DIR = bin
SRC_DIR = src
OBJ_DIR = build
MOD_DIR = include
LIB_DIR = lib

# List of source files
SRCS = $(wildcard $(SRC_DIR)/*.f90)
LIBS = ${MKLROOT}/include/mkl_spblas.f90 lib/Quantum_Ising_1D.f90
# $(wildcard $(LIB_DIR)/*.f90)

# Generate object file names based on source file names
SRC_OBJS = $(patsubst $(SRC_DIR)/%.f90,$(EXE_DIR)/%,$(SRCS))
# Separate LIBS into two parts based on association
LIBS_LIB_DIR = $(filter $(LIB_DIR)/%,$(LIBS))
LIBS_MKL_DIR = $(filter ${MKLROOT}/include/%,$(LIBS))

# Generate object file names based on library file names
LIB_OBJS = $(patsubst ${MKLROOT}/include/%.f90,$(OBJ_DIR)/%.o,$(LIBS_MKL_DIR))\
    $(patsubst $(LIB_DIR)/%.f90,$(OBJ_DIR)/%.o,$(LIBS_LIB_DIR))

# Target to compile all object files and generate module files
all: $(LIB_OBJS) $(SRC_OBJS)
	
# Rule to compile individual source files into object files
$(EXE_DIR)/%: $(SRC_DIR)/%.f90 $(LIB_OBJS)
	$(FC) $(FFLAGS) -module include -o $@ -I include ${MKL_INCLUDE} $< $(LIB_OBJS) ${STATICLIB}

# Rule to compile individual source files into object files
$(OBJ_DIR)/%.o: $(LIB_DIR)/%.f90
	$(FC) $(FFLAGS) -c -module include -o $@ $<

$(OBJ_DIR)/%.o: ${MKLROOT}/include/%.f90
	$(FC) $(FFLAGS) -c -module include -o $@ $<

clean:
	rm -f build/* include/* bin/*

.PHONY: clean
