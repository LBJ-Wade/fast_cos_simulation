///
/// \file  comm.c
/// \brief Functions for MPI communications
///

#include <assert.h>

#ifdef MPI
#include <mpi.h>
#else
#include <stdlib.h>
#endif

static int this_node= -1;
static int n_nodes= 0;
static int parallel_level= 0;
  // 0: no MPI
  // 1: MPI_THREAD_SINGLE, only one thread per MPI node
  // 2: MPI_THREAD_FUNNELED, only the thread that called MPI_Init_thread will
  //    make MPI calls.

#ifdef MPI

///
/// Initialisation
/// abc
///

void comm_mpi_init(int* p_argc, char*** p_argv)
{
#ifdef _OPENMP
  int thread_level;
  MPI_Init_thread(p_argc, p_argv, MPI_THREAD_FUNNELED, &thread_level);
  if(thread_level >= MPI_THREAD_FUNNELED)
    parallel_level= 2;
  else
    parallel_level= 1;
#else
  int ret= MPI_Init(p_argc, p_argv); assert(ret == MPI_SUCCESS);
  parallel_level= 1;
#endif
  

  MPI_Comm_rank(MPI_COMM_WORLD, &this_node);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
}

void comm_mpi_finalise(void)
{
  MPI_Finalize();
}

void comm_abort(void)
{
  MPI_Abort(MPI_COMM_WORLD, 1);
}


void comm_bcast_int(int* p_int, int count)
{
  MPI_Bcast(p_int, count, MPI_INT, 0, MPI_COMM_WORLD);
}

void comm_bcast_double(double* p_double, int count)
{
  MPI_Bcast(p_double, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

#else

//
// Serial version
//

void comm_mpi_init(int* p_argc, char*** p_argv)
{
  this_node= 0;
  parallel_level= 0;
}

void comm_mpi_finalise(void)
{
}

void comm_abort(void)
{
  abort();
}


void comm_bcast_int(int* p_int, int count)
{
}

void comm_bcast_double(double* p_double, int count)
{
}

#endif

//
// Function usable for both serial and MPI version
//

///
/// Initialisation
/// abc
///

int comm_this_node(void)
{
  return this_node;
}

int comm_n_nodes(void)
{
  return n_nodes;
}
