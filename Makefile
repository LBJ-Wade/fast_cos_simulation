#
# fs: fast cosmological simulation system
#   

#
# Compile options
#

#OPT+= -DDOUBLEPRECISION
OPT += -DMPI

#
# Compile configurations
#

# Define OPENMP to enable MPI+OpenMP hybrid parallelization
# OPENMP  = -fopenmp # -openmp for Intel, -fopenmp for gcc, llvm doesn't support

ifeq (,$(findstring -DMPI, $(OPT)))
CC       = cc
else
CC       = mpicc -std=c99
endif

WOPT    ?= -Wall
CFLAGS  := -O3 $(WOPT) $(OPENMP) $(OPT)
LIBS    := -lm

# Define paths of FFTW3 & GSL libraries if necessary.

LUA_DIR   ?= #e.g. /opt/local
FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH = $(FFTW3_DIR) $(GSL_DIR) $(LUA_DIR)

CFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS   += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

EXEC = fs libfs.a # doc
all: $(EXEC)

OBJS := main.o comm.o msg.o power.o cosmology.o mem.o util.o fft.o config.o
OBJS += lpt.o pm.o cola.o

#
# Linking libraries
#
LIBS += -llua -ldl 
LIBS += -lgsl -lgslcblas

ifeq (,$(findstring -DDOUBLEPRECISION, $(OPT)))
  # Single precision FFTW
  FFTWSUF=f
endif
LIBS += -lfftw3$(FFTWSUF)

ifeq (-DMPI,$(findstring -DMPI, $(OPT)))
  LIBS += -lfftw3$(FFTWSUF)_mpi
endif

ifdef OPENMP
  LIBS += -lfftw3$(FFTWSUF)_omp
  #LIBS += -lfftw3$(FFTWSUF)_threads # for thread parallelization instead of omp
endif

# Compiling rule
fs: $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $@

# Library libfs.a
libfs.a: $(OBJS)
	ar r $@ $(OBJS)

doc: $(OBJS)
	cd doc && doxygen >& doxygen.log


.PHONY: clean run dependence
clean:
	rm -f $(EXEC) $(OBJS)

run:
	mpirun -n 2 fs

dependence:
	gcc -MM -MG *.c

