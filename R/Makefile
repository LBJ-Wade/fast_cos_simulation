#
# Rfs: R library
#   
.SECONDARY:

#
# Compile options
#

#OPT+= -DDOUBLEPRECISION
# R library is non-MPI

#
# Compile configurations
#

# Define OPENMP to enable OpenMP parallelization
# OPENMP  = -fopenmp # -openmp for Intel, -fopenmp for gcc, llvm doesn't support

CC       = cc

WOPT    ?= -Wall
CFLAGS  := -O3 $(WOPT) $(OPENMP) $(OPT) -I../
LIBS    := -lm

# Define paths of FFTW3 & GSL libraries if necessary.

LUA_DIR   ?= #e.g. /opt/local
FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH = $(FFTW3_DIR) $(GSL_DIR) $(LUA_DIR)

CFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS   += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

EXEC = Rfs.so
all: $(EXEC)

OBJS += comm.o msg.o power.o cosmology.o mem.o util.o fft.o config.o
OBJS += lpt.o


#
# Linking libraries
#
#LIBS += -llua -ldl 
LIBS += -lgsl -lgslcblas

ifeq (,$(findstring -DDOUBLEPRECISION, $(OPT)))
  FFTWSUF=f    # Single precision FFTW
endif
LIBS += -lfftw3$(FFTWSUF)

ifdef OPENMP
  LIBS += -lfftw3$(FFTWSUF)_omp
  #LIBS += -lfftw3$(FFTWSUF)_threads # for thread parallelization instead of omp
endif


# Compiling rule
Rfs.so: Rfs.c $(OBJS)
	MAKEFLAGS="CFLAGS=-I../" R CMD SHLIB Rfs.c $(OBJS) $(LIBS)

%.c: ../%.c
	ln -s $< $@

.PHONY: clean run dependence
clean:
	rm -f $(EXEC) $(OBJS) Rfs.o

dependence:
	gcc -MM -MG *.c
