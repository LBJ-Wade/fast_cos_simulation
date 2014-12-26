#
# fs
#   

# Define OPENMP to enable MPI+OpenMP hybrid parallelization
# OPENMP  = -fopenmp # -openmp for Intel, -fopenmp for gcc

CC      = mpicc -std=c99 
WOPT    ?= -Wall
CFLAGS  := -O3 $(WOPT) $(OPENMP)
LIBS    := -lm


# Compile options


# Define paths of FFTW3 & GSL libraries if necessary.

LUA_DIR   ?= #e.g. /opt/local
FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH = $(FFTW3_DIR) $(GSL_DIR) $(LUA_DIR)

CFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS   += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

EXEC = fs
all: $(EXEC)

OBJS := main.o comm.o msg.o

LIBS += -llua -ldl 
LIBS += -lgsl -lgslcblas
LIBS += -lfftw3f_mpi -lfftw3f


ifdef OPENMP
  LIBS += -lfftw3f_omp
  #LIBS += -lfftw3f_threads       # for thread parallelization instead of omp
endif

fs: $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $@

.PHONY: clean run dependence
clean:
	rm -f $(EXEC) $(OBJS) $(OBJS2) move_min.?

run:
	mpirun -n 2 fs

dependence:
	gcc -MM -MG *.c
