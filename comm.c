//
// Taking care of MPI communications
//

#include <assert.h>
#include <mpi.h>

static int this_node= -1;

int comm_mpi_init(int* p_argc, char*** p_argv)
{
  int ret= MPI_Init(p_argc, p_argv); assert(ret == MPI_SUCCESS);
  MPI_Comm_rank(MPI_COMM_WORLD, &this_node);
  
  return 0;
}

int comm_this_node(void)
{
  return this_node;
}
